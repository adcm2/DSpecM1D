#pragma once

#include "pole_response_spool.hpp"
#include <GSHTrans/Utility.hpp>
#include <GaussQuad/All>

#include <chrono>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <numeric>
#include <set>

namespace dspec_global_gsh {
struct ConsumerOptions {
  std::size_t sampleAngleCount=0; // zero means the complete bandwidth grid
  std::size_t longitudeTile=128;
  std::size_t workers=1;
  std::size_t maxRadialSides=8;
  std::vector<std::size_t> radialSideIndices;
  double maxOutputGiB=8.0;
  double maxWorkingGiB=8.0;
  double transformScratchGiBPerWorker=2.0;
};
struct ConsumerResult {
  std::size_t angularPointCount=0,fullAngularPointCount=0,frequencyCount=0,radialSideCount=0;
  std::uint64_t outputBytes=0;
  double elapsedSeconds=0.0;
  double serialContextSetupSeconds=0.0,workerExecutionWallSeconds=0.0;
};

/** Replays compact source/unit-adjoint pole responses through full-band GSH
 * inverse transforms. A native-side accumulator completes all frequency
 * products before output chunks are written once. */
inline ConsumerResult runPoleResponseSpoolKernel(
    const std::string &spoolPath,const Full1D::SEM &sem,
    const std::string &outputPath,ConsumerOptions options={}) {
  using Clock=std::chrono::steady_clock;
  const auto functionStarted=Clock::now();
  if(options.workers==0||options.workers>48||options.longitudeTile==0||
     !(std::isfinite(options.maxOutputGiB)&&options.maxOutputGiB>0.0&&
       std::isfinite(options.maxWorkingGiB)&&options.maxWorkingGiB>0.0&&
       std::isfinite(options.transformScratchGiBPerWorker)&&options.transformScratchGiBPerWorker>0.0))
    throw std::invalid_argument("invalid global kernel consumer worker, tile, or memory limits");
  PoleResponseSpoolReader metadata(spoolPath);
  const int lmax=metadata.lmax();
  if(lmax>1536)throw std::runtime_error("global kernel consumer exceeds the angular-library degree ceiling");
  const std::size_t nTheta=static_cast<std::size_t>(2*lmax+1);
  const std::size_t nPhi=static_cast<std::size_t>(GSHTrans::FastFFTSize(4*lmax+1));
  if(nTheta>std::numeric_limits<std::size_t>::max()/nPhi)
    throw std::overflow_error("global angular point count overflows size_t");
  const std::size_t fullAngles=nTheta*nPhi;
  const std::size_t sampleAngles=options.sampleAngleCount==0?fullAngles:options.sampleAngleCount;
  if(sampleAngles==0||sampleAngles>fullAngles)
    throw std::invalid_argument("requested angular prefix is outside the complete bandwidth grid");
  const auto &sides=metadata.sides();
  if(sides.empty()||options.maxRadialSides==0)
    throw std::invalid_argument("global kernel consumer requires native sides and a positive side cap");
  std::vector<std::size_t> selectedSides=options.radialSideIndices;
  if(selectedSides.empty()) {
    if(sides.size()>options.maxRadialSides)
      throw std::invalid_argument("explicit native-side indices are required when the spool exceeds the side cap");
    selectedSides.resize(sides.size());std::iota(selectedSides.begin(),selectedSides.end(),std::size_t{0});
  }
  if(selectedSides.size()>options.maxRadialSides||selectedSides.empty())
    throw std::invalid_argument("selected native-side count is outside the configured cap");
  std::set<std::size_t> uniqueSides;
  for(const auto index:selectedSides)
    if(index>=sides.size()||!uniqueSides.insert(index).second)
      throw std::invalid_argument("selected native-side indices are invalid or duplicated");
  if(sampleAngles>std::numeric_limits<std::size_t>::max()/selectedSides.size()/6/sizeof(double))
    throw std::overflow_error("global kernel output size overflows size_t");
  const std::uint64_t bytes=static_cast<std::uint64_t>(sampleAngles)*selectedSides.size()*6*sizeof(double);
  constexpr double gib=1024.0*1024.0*1024.0;
  if(options.maxOutputGiB>static_cast<double>(std::numeric_limits<std::uint64_t>::max())/gib||
     options.maxWorkingGiB>static_cast<double>(std::numeric_limits<std::uint64_t>::max())/gib||
     options.transformScratchGiBPerWorker>static_cast<double>(std::numeric_limits<std::uint64_t>::max())/gib||
     bytes>static_cast<std::uint64_t>(options.maxOutputGiB*gib))
    throw std::runtime_error("global kernel output exceeds configured byte cap");
  const auto parent=std::filesystem::path(outputPath).parent_path();
  if(!parent.empty())std::filesystem::create_directories(parent);
  if(std::filesystem::exists(outputPath)||std::filesystem::exists(outputPath+".json"))
    throw std::runtime_error("global kernel output or manifest path already exists");
  const auto space=std::filesystem::space(parent.empty()?std::filesystem::current_path():parent);
  if(bytes>space.available)throw std::runtime_error("insufficient free disk space for global kernel output");

  const std::size_t fullSideAccumulatorBytes=fullAngles*6*sizeof(double);
  const std::size_t rotationBytes=2*SPARSESPEC::Windowed3DFieldEngine::angularWignerBytesPerAngle(lmax,lmax);
  const std::size_t generatedAngularScratchBytes=SPARSESPEC::Windowed3DFieldEngine::angularTableBytesPerAngle(lmax,lmax);
  const std::size_t readerPersistentBytes=sides.size()*(2*sizeof(std::int32_t)+3*sizeof(double)+sizeof(std::uint8_t))+
      metadata.frequencyCount()*(sizeof(std::int32_t)+2*sizeof(Complex)+sizeof(double));
  const std::size_t coefficientCount=12*static_cast<std::size_t>(lmax+1)*static_cast<std::size_t>(lmax+1);
  const std::size_t coefficientBytes=coefficientCount*sizeof(Complex);
  const std::uint64_t transformScratchBytes=static_cast<std::uint64_t>(options.transformScratchGiBPerWorker*gib);
  const std::size_t lowDegreeBatchScratch=lmax<=32?
      4*(fullAngles+static_cast<std::size_t>(lmax+1)*static_cast<std::size_t>(lmax+1))*sizeof(Complex):0;
  const std::size_t knownPerWorkerBytes=fullAngles*(12*sizeof(Complex)+6*sizeof(double))+fullSideAccumulatorBytes+
      rotationBytes+generatedAngularScratchBytes+coefficientBytes+lowDegreeBatchScratch+readerPersistentBytes;
  if(transformScratchBytes>std::numeric_limits<std::size_t>::max()-knownPerWorkerBytes)
    throw std::overflow_error("global transform scratch reserve overflows addressable memory");
  const std::size_t perWorkerBytes=knownPerWorkerBytes+static_cast<std::size_t>(transformScratchBytes);
  const std::uint64_t maxWorkingBytes=static_cast<std::uint64_t>(options.maxWorkingGiB*gib);
  if(sides.size()>std::numeric_limits<std::size_t>::max()/metadata.frequencyCount()/4/16)
    throw std::overflow_error("spool offset validation working memory overflows");
  const std::size_t spoolValidationBytes=metadata.frequencyCount()*sides.size()*4*16;
  if(spoolValidationBytes>256u*1024u*1024u||maxWorkingBytes<=spoolValidationBytes)
    throw std::runtime_error("spool validation memory exceeds its cap or the configured working-memory budget");
  const std::size_t workersByMemory=perWorkerBytes==0?0:
      static_cast<std::size_t>((maxWorkingBytes-spoolValidationBytes)/perWorkerBytes);
  const std::size_t workers=std::min({options.workers,selectedSides.size(),static_cast<std::size_t>(48),workersByMemory});
  if(workers==0)throw std::runtime_error("configured working-memory cap cannot hold one complete-band side/frequency workspace");
  {
    std::ofstream initial(outputPath,std::ios::binary|std::ios::trunc);
    if(!initial)throw std::runtime_error("cannot create global kernel output");
    initial.seekp(static_cast<std::streamoff>(bytes-1));initial.put('\0');initial.close();
    if(!initial)throw std::runtime_error("cannot size global kernel output");
  }

  struct Worker {
    std::unique_ptr<PoleResponseSpoolReader> reader;
    std::unique_ptr<GlobalGridPlan> plan;
    std::fstream output;
    double sideSeconds=0.0,recordReadSeconds=0.0,frequencyAccumulationSeconds=0.0,outputWriteSeconds=0.0;
  };
  std::vector<Worker> contexts(workers);
  Eigen::setNbThreads(1);
  const auto setupStarted=Clock::now();
  for(auto &ctx:contexts) {
    ctx.reader=std::make_unique<PoleResponseSpoolReader>(spoolPath);
    ctx.plan=std::make_unique<GlobalGridPlan>(lmax);
    ctx.output.open(outputPath,std::ios::binary|std::ios::in|std::ios::out);
    if(!ctx.output)throw std::runtime_error("worker cannot open global kernel output");
  }
  const double setupSeconds=std::chrono::duration<double>(Clock::now()-setupStarted).count();
  const auto workerWallStarted=Clock::now();
  runOwnedTiles(selectedSides.size(),workers,[&](std::size_t outputSide,std::size_t worker) {
    auto &ctx=contexts[worker];
    const auto workStart=std::chrono::steady_clock::now();
    std::vector<std::array<double,6>> accumulator(fullAngles);
    for(std::size_t bin=0;bin<ctx.reader->frequencyCount();++bin) {
      const auto readStarted=Clock::now();
      auto records=ctx.reader->read(bin,selectedSides[outputSide]);
      ctx.recordReadSeconds+=std::chrono::duration<double>(Clock::now()-readStarted).count();
      if(records.audit.basisReconstructionRelative>1e-12)
        throw std::runtime_error("receiver unit-basis reconstruction exceeds the 1e-12 gate");
      if(records.audit.adjointMixing==Complex{})continue;
      const auto frequencyKernel=ctx.plan->evaluate(records.forward,records.unitAdjoint,
          records.l0,records.audit.adjointMixing,sem);
      if(frequencyKernel.size()!=fullAngles)throw std::logic_error("complete-band transform returned the wrong grid size");
      const auto accumulationStarted=Clock::now();
      for(std::size_t p=0;p<fullAngles;++p)for(std::size_t c=0;c<6;++c)accumulator[p][c]+=frequencyKernel[p][c];
      ctx.frequencyAccumulationSeconds+=std::chrono::duration<double>(Clock::now()-accumulationStarted).count();
    }
    const auto writeStarted=Clock::now();
    for(std::size_t begin=0;begin<sampleAngles;begin+=options.longitudeTile) {
      const std::size_t count=std::min(options.longitudeTile,sampleAngles-begin);
      const std::uint64_t offset=(static_cast<std::uint64_t>(outputSide)*sampleAngles+begin)*6*sizeof(double);
      ctx.output.clear();ctx.output.seekp(static_cast<std::streamoff>(offset));
      ctx.output.write(reinterpret_cast<const char*>(accumulator.data()+begin),static_cast<std::streamsize>(count*sizeof(accumulator.front())));
      if(!ctx.output)throw std::runtime_error("global kernel side-tile write failed");
    }
    ctx.outputWriteSeconds+=std::chrono::duration<double>(Clock::now()-writeStarted).count();
    ctx.sideSeconds+=std::chrono::duration<double>(Clock::now()-workStart).count();
  });

  const double workerWallSeconds=std::chrono::duration<double>(Clock::now()-workerWallStarted).count();
  for(auto &ctx:contexts)if(ctx.output.is_open())ctx.output.close();
  const double elapsed=std::chrono::duration<double>(Clock::now()-functionStarted).count();
  double workerElapsedSum=0.0,readSum=0.0,rotationSum=0.0,inverseSum=0.0,contractionSum=0.0,accumulationSum=0.0,writeSum=0.0;
  for(const auto &ctx:contexts){
    workerElapsedSum+=ctx.sideSeconds;readSum+=ctx.recordReadSeconds;
    rotationSum+=ctx.plan->timings().coefficientRotationLayoutSeconds;
    inverseSum+=ctx.plan->timings().inverseGshFftSeconds;
    contractionSum+=ctx.plan->timings().pointKernelSeconds;
    accumulationSum+=ctx.frequencyAccumulationSeconds;writeSum+=ctx.outputWriteSeconds;
  }
  std::ofstream manifest(outputPath+".json");manifest<<std::setprecision(17);
  manifest<<"{\"mode\":\"compact_response_global_gsh_fft_frequency_accumulator\",\"degree\":"<<lmax
    <<",\"angular_rule\":\"ForBand(L,2,2.0) nodes and FastFFTSize(4L+1) longitudes; band-L global GSH synthesis\""
    <<",\"n_theta\":"<<nTheta<<",\"n_phi\":"<<nPhi<<",\"full_angular_points\":"<<fullAngles
    <<",\"sampled_angular_points\":"<<sampleAngles<<",\"angular_subset_is_complete\":"<<(sampleAngles==fullAngles?"true":"false")
    <<",\"radial_side_count\":"<<selectedSides.size()<<",\"native_side_indices\":[";
  for(std::size_t i=0;i<selectedSides.size();++i){if(i)manifest<<',';manifest<<selectedSides[i];}
  manifest<<"],\"native_side_metadata\":[";
  for(std::size_t i=0;i<selectedSides.size();++i){
    if(i)manifest<<',';const auto &side=sides[selectedSides[i]];
    manifest<<"{\"element\":"<<side.element<<",\"node\":"<<side.node<<",\"radius_m\":"<<side.radiusMetres
      <<",\"radius_model\":"<<side.radiusModel<<",\"width_m\":"<<side.widthMetres<<",\"fluid\":"<<(side.fluid?"true":"false")<<'}';
  }
  manifest<<"],\"frequency_count\":"<<metadata.frequencyCount()
    <<",\"frequency_records\":[";
  for(std::size_t bin=0;bin<metadata.frequencyCount();++bin){
    if(bin)manifest<<',';
    const auto &a=metadata.frequencyAudit(bin);
    manifest<<"{\"spool_index\":"<<bin<<",\"original_bin\":"<<a.frequencyBin
      <<",\"z\":["<<a.z.real()<<','<<a.z.imag()<<"],\"adjoint_mixing\":["
      <<a.adjointMixing.real()<<','<<a.adjointMixing.imag()<<"],\"receiver_basis_reconstruction_relative\":"
      <<a.basisReconstructionRelative<<'}';
  }
  manifest<<"]"
    <<",\"angular_transform_backend\":\"GSHTrans::TransformKernel::Loop with FFT longitude synthesis; same-spin channel batches at L<=32; Eigen GEMM is not used for angular synthesis\""
    <<",\"output_write_chunk_points\":"<<options.longitudeTile<<",\"workers\":"<<workers
    <<",\"max_radial_sides_cap\":"<<options.maxRadialSides<<",\"configured_working_memory_gib\":"<<options.maxWorkingGiB
    <<",\"known_working_bytes_per_worker\":"<<knownPerWorkerBytes
    <<",\"reader_persistent_metadata_bytes_per_worker\":"<<readerPersistentBytes
    <<",\"configured_transform_scratch_reserve_gib_per_worker\":"<<options.transformScratchGiBPerWorker
    <<",\"estimated_working_bytes_per_worker\":"<<perWorkerBytes
    <<",\"spool_offset_validation_scratch_bytes\":"<<spoolValidationBytes
    <<",\"output_bytes\":"<<bytes<<",\"output_layout\":\"radial-side-major; theta-major angular prefix; point-major channels A,C,F,L,N,kappa; complete matching-frequency products accumulated per native side before output\""
    <<",\"kernel_is_full_earth_volume_integral\":false,\"source_contacts_included\":false"
    <<",\"binary_file\":\""<<std::filesystem::path(outputPath).filename().string()<<"\",\"elapsed_s\":"<<elapsed
    <<",\"timings\":{\"serial_context_setup_wall_s\":"<<setupSeconds<<",\"worker_execution_wall_s\":"<<workerWallSeconds
    <<",\"total_wall_s\":"<<elapsed<<",\"worker_elapsed_sum_s\":"<<workerElapsedSum
    <<",\"worker_stage_elapsed_sums\":{\"spool_reads\":"<<readSum<<",\"coefficient_rotation_layout\":"<<rotationSum
    <<",\"gsh_inverse_fft\":"<<inverseSum<<",\"point_kernel_contractions\":"<<contractionSum
    <<",\"frequency_accumulation\":"<<accumulationSum<<",\"output_writes\":"<<writeSum<<"},\"workers\":[";
  for(std::size_t i=0;i<contexts.size();++i){if(i)manifest<<',';const auto &ctx=contexts[i];const auto &t=ctx.plan->timings();
    manifest<<"{\"side_elapsed_s\":"<<ctx.sideSeconds<<",\"spool_reads_s\":"<<ctx.recordReadSeconds
      <<",\"coefficient_rotation_layout_s\":"<<t.coefficientRotationLayoutSeconds
      <<",\"gsh_inverse_fft_s\":"<<t.inverseGshFftSeconds<<",\"point_kernel_contractions_s\":"<<t.pointKernelSeconds
      <<",\"frequency_accumulation_s\":"<<ctx.frequencyAccumulationSeconds<<",\"output_writes_s\":"<<ctx.outputWriteSeconds<<'}';}
  manifest<<"],\"worker_timing_note\":\"Per-worker durations overlap in wall time and their sums are not a non-overlapping CPU-time measure. Total wall includes serial context setup and worker execution.\""
    <<",\"worker_owned_side_seconds\":[";
  for(std::size_t i=0;i<contexts.size();++i){if(i)manifest<<',';manifest<<contexts[i].sideSeconds;}
  manifest<<"],\"known_side_accumulator_bytes\":"<<fullSideAccumulatorBytes
    <<",\"memory_note\":\"Preflight counts twelve complex fields, a frequency-result buffer, the full-side accumulator, twelve coefficient arrays, two frame Wigners, one Generated-mode angular scratch table, and low-degree batch scratch where applicable, plus a configured transform scratch reserve. Generated-mode Wigner work is treated as chunked rather than retained for every theta row. Actual FFTW/GSHTrans/Eigen/SEM/filesystem/thread-stack peak use must be measured; this is not a peak-RSS guarantee.\"}\n";
  if(!manifest)throw std::runtime_error("failed writing global kernel manifest");
  return {sampleAngles,fullAngles,metadata.frequencyCount(),selectedSides.size(),bytes,elapsed,setupSeconds,workerWallSeconds};
}
} // namespace dspec_global_gsh
