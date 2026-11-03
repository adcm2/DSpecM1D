#pragma once

#include "pole_response_spool.hpp"
#include <GSHTrans/Utility.hpp>
#include <GSHTrans/OpenMP.hpp>
#include <GaussQuad/All>
#ifdef _OPENMP
#include <omp.h>
#endif

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
namespace consumer_detail {
inline std::size_t checkedMultiply(std::size_t a,std::size_t b,const char *what) {
  if(b!=0&&a>std::numeric_limits<std::size_t>::max()/b)throw std::overflow_error(what);
  return a*b;
}
inline std::size_t checkedAdd(std::size_t a,std::size_t b,const char *what) {
  if(a>std::numeric_limits<std::size_t>::max()-b)throw std::overflow_error(what);
  return a+b;
}
inline std::size_t angularWignerBytesPerAngle(int maxDegree,int mMax) {
  if(maxDegree<0)throw std::invalid_argument("angular memory estimate requires a nonnegative degree");
  const int nMax=std::min(2,maxDegree);mMax=std::min(mMax,maxDegree);
  std::size_t values=0;
  for(int n=-nMax;n<=nMax;++n)
    values=checkedAdd(values,static_cast<std::size_t>(GSHTrans::GSHIndices<GSHTrans::All>(maxDegree,mMax,n).Size()),
        "angular Wigner table size overflows");
  const auto payload=checkedMultiply(values,sizeof(double),"angular Wigner table bytes overflow");
  const auto indices=checkedMultiply(static_cast<std::size_t>(2*nMax+1),sizeof(std::size_t),
      "angular Wigner index bytes overflow");
  return checkedAdd(payload,indices,"angular Wigner table total bytes overflow");
}
inline std::size_t angularTableBytesPerAngle(int maxDegree,int mMax) {
  if(maxDegree<0)throw std::invalid_argument("angular memory estimate requires a nonnegative degree");
  const std::size_t orders=static_cast<std::size_t>(2*std::min(2,maxDegree)+1);
  const auto coefficients=checkedMultiply(checkedMultiply(orders,5u,"angular coefficient table size overflows"),
      static_cast<std::size_t>(maxDegree+1),"angular coefficient table size overflows");
  const auto coefficientBytes=checkedMultiply(coefficients,sizeof(double),"angular coefficient table bytes overflow");
  return checkedAdd(angularWignerBytesPerAngle(maxDegree,mMax),coefficientBytes,"angular table total bytes overflow");
}
} // namespace consumer_detail

struct ConsumerOptions {
  std::size_t sampleAngleCount=0; // zero means the complete bandwidth grid
  std::size_t longitudeTile=128;
  std::size_t workers=1;
  std::size_t gshThreads=1;
  std::size_t maxRadialSides=8;
  bool storedWignerValues=false;
  std::vector<std::size_t> radialSideIndices;
  double maxOutputGiB=8.0;
  double maxWorkingGiB=8.0;
  double transformScratchGiBPerWorker=2.0;
};
struct ConsumerResult {
  std::size_t angularPointCount=0,fullAngularPointCount=0,frequencyCount=0,radialSideCount=0;
  std::uint64_t outputBytes=0;
  double elapsedSeconds=0.0;
  double serialGridContextSetupSeconds=0.0,workerExecutionWallSeconds=0.0;
  std::uint64_t sharedWignerBytes=0,estimatedWorkingBytes=0;
  std::uint64_t frequencySideTransformsExecuted=0,frequencySideZeroMixingSkipped=0;
  std::size_t sameSpinBatchSize=1;
  std::size_t gshThreads=1,workerCount=0,requestedThreadProduct=0,effectiveThreadProductBound=0;
  std::uint64_t innerThreadScratchBytesPerWorker=0,innerThreadScratchBytesTotal=0;
};

/** Exact double-table extent for this grid: (2L+1) theta rows times the
 * sum of GSHIndices<All>(L,L,n).Size() for n=-2..2. */
inline std::uint64_t storedWignerTableBytes(int lmax) {
  if(lmax<1||lmax>1536)throw std::invalid_argument("stored Wigner table degree is outside the supported range");
  const std::uint64_t side=static_cast<std::uint64_t>(lmax)+1;
  const std::uint64_t coefficients=5*side*side-10;
  return static_cast<std::uint64_t>(2*lmax+1)*coefficients*sizeof(double);
}

/** Replays compact source/unit-adjoint pole responses through full-band GSH
 * inverse transforms. A native-side accumulator completes all frequency
 * products before output chunks are written once. */
inline ConsumerResult runPoleResponseSpoolKernel(
    const std::string &spoolPath,const Full1D::SEM &sem,
    const std::string &outputPath,ConsumerOptions options={}) {
  using Clock=std::chrono::steady_clock;
  const auto functionStarted=Clock::now();
  if(options.workers==0||options.workers>48||options.gshThreads==0||options.gshThreads>48||
     options.workers>48/options.gshThreads||options.longitudeTile==0||
     !(std::isfinite(options.maxOutputGiB)&&options.maxOutputGiB>0.0&&
       std::isfinite(options.maxWorkingGiB)&&options.maxWorkingGiB>0.0&&
       std::isfinite(options.transformScratchGiBPerWorker)&&options.transformScratchGiBPerWorker>0.0))
    throw std::invalid_argument("invalid global kernel consumer worker, tile, or memory limits");
  if(options.gshThreads>1&&!GSHTrans::OpenMP::Available)
    throw std::runtime_error("--global-gsh-threads greater than one requires an OpenMP-enabled GSHTrans build");
#ifdef _OPENMP
  if(options.gshThreads>1)omp_set_dynamic(0);
#endif
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
  const std::size_t rotationBytes=consumer_detail::checkedMultiply(2,consumer_detail::angularWignerBytesPerAngle(lmax,lmax),
      "per-worker frame rotation tables overflow");
  const std::size_t framePhaseCacheBytes=2*(static_cast<std::size_t>(lmax+1)*static_cast<std::size_t>(lmax+1)*sizeof(Complex)+
      static_cast<std::size_t>(lmax+1)*sizeof(double));
  const std::size_t generatedAngularScratchBytes=options.storedWignerValues?0:
      consumer_detail::angularTableBytesPerAngle(lmax,lmax);
  const std::uint64_t sharedWignerBytes=options.storedWignerValues?storedWignerTableBytes(lmax):0;
  const std::size_t sideVectorBytes=consumer_detail::checkedMultiply(
      consumer_detail::checkedAdd(sides.size(),1,"reader side-vector capacity overflows"),sizeof(NativeSide),
      "reader side-vector bytes overflow");
  const std::size_t auditVectorBytes=consumer_detail::checkedMultiply(
      consumer_detail::checkedAdd(metadata.frequencyCount(),1,"reader audit-vector capacity overflows"),sizeof(SpoolAudit),
      "reader audit-vector bytes overflow");
  const std::size_t readerPersistentBytes=consumer_detail::checkedAdd(sizeof(PoleResponseSpoolReader),
      consumer_detail::checkedAdd(sideVectorBytes,auditVectorBytes,"reader metadata bytes overflow"),
      "reader instance bytes overflow");
  // At degree l, recorded local orders are -min(2,l)..min(2,l): one
  // coefficient at l=0, three at l=1, then five per degree, hence 5L-1.
  const std::size_t maxPoleLabels=static_cast<std::size_t>(5*lmax-1);
  const std::size_t decodedLabelBytes=consumer_detail::checkedAdd(sizeof(std::pair<int,int>),
      sizeof(std::array<Complex,6>),"decoded pole label size overflows");
  const std::size_t oneDecodedVectorBytes=consumer_detail::checkedMultiply(
      consumer_detail::checkedAdd(maxPoleLabels,1,"decoded pole vector capacity overflows"),decodedLabelBytes,
      "decoded pole vector bytes overflow");
  const std::size_t decodedRecordsBytes=consumer_detail::checkedAdd(sizeof(FrequencySideRecords),
      consumer_detail::checkedMultiply(2,oneDecodedVectorBytes,"decoded pole pair bytes overflow"),
      "decoded frequency-side record bytes overflow");
  const std::size_t coefficientCount=12*static_cast<std::size_t>(lmax+1)*static_cast<std::size_t>(lmax+1);
  const std::size_t coefficientBytes=coefficientCount*sizeof(Complex);
  const std::size_t maxTransformBatch=lmax<=32?4:2;
  // Thread-local workspaces are cached by batch count. Low-degree plans use
  // both four-channel spin-0 batches and two-channel other-spin batches;
  // high-degree plans may retain both their selected batch and a one-channel
  // tail shape, so charge both shapes conservatively.
  const std::size_t threadLocalBatchStorage=maxTransformBatch+(maxTransformBatch==4?2:1);
  const std::size_t generatedWignerTlsBytes=options.storedWignerValues?0:
      consumer_detail::checkedMultiply(consumer_detail::checkedMultiply(static_cast<std::size_t>(lmax+1),
          static_cast<std::size_t>(lmax+1),"Generated Wigner TLS extent overflows"),sizeof(double),
          "Generated Wigner TLS bytes overflow");
  const std::size_t fftBufferComplexCount=consumer_detail::checkedMultiply(
      consumer_detail::checkedMultiply(2,threadLocalBatchStorage,"per-inner-thread FFT buffer extent overflows"),
      nPhi,"per-inner-thread FFT buffer extent overflows");
  const std::size_t fftBuffersTlsBytes=consumer_detail::checkedMultiply(fftBufferComplexCount,sizeof(Complex),
      "per-inner-thread FFT buffer bytes overflow");
  const std::size_t knownTlsBytesPerInnerThread=consumer_detail::checkedAdd(generatedWignerTlsBytes,fftBuffersTlsBytes,
      "per-inner-thread known TLS bytes overflow");
  const std::size_t innerThreadScratchBytesPerWorker=consumer_detail::checkedMultiply(options.gshThreads,
      knownTlsBytesPerInnerThread,"per-worker inner-thread TLS bytes overflow");
  const std::uint64_t transformScratchBytes=static_cast<std::uint64_t>(options.transformScratchGiBPerWorker*gib);
  const std::size_t basePerWorkerBytes=fullAngles*(12*sizeof(Complex)+6*sizeof(double))+fullSideAccumulatorBytes+
      rotationBytes+framePhaseCacheBytes+generatedAngularScratchBytes+coefficientBytes+readerPersistentBytes+decodedRecordsBytes+
      innerThreadScratchBytesPerWorker;
  const std::size_t batchTwoScratch=2*(fullAngles+static_cast<std::size_t>(lmax+1)*static_cast<std::size_t>(lmax+1))*sizeof(Complex);
  const std::size_t batchFourScratch=4*(fullAngles+static_cast<std::size_t>(lmax+1)*static_cast<std::size_t>(lmax+1))*sizeof(Complex);
  if(transformScratchBytes>std::numeric_limits<std::size_t>::max()-basePerWorkerBytes)
    throw std::overflow_error("global transform scratch reserve overflows addressable memory");
  const std::size_t basePerWorkerWithReserve=basePerWorkerBytes+static_cast<std::size_t>(transformScratchBytes);
  const std::uint64_t maxWorkingBytes=static_cast<std::uint64_t>(options.maxWorkingGiB*gib);
  if(sides.size()>std::numeric_limits<std::size_t>::max()/metadata.frequencyCount()/4/16)
    throw std::overflow_error("spool offset validation working memory overflows");
  const std::size_t spoolValidationBytes=consumer_detail::checkedMultiply(
      consumer_detail::checkedMultiply(metadata.frequencyCount(),sides.size(),"spool validation key count overflows"),
      4*sizeof(std::uint64_t)*2,"spool validation bytes overflow");
  if(spoolValidationBytes>256u*1024u*1024u||maxWorkingBytes<=spoolValidationBytes||
     sharedWignerBytes>maxWorkingBytes-spoolValidationBytes)
    throw std::runtime_error("spool validation memory exceeds its cap or the configured working-memory budget");
  const std::uint64_t sharedResidentBytes=consumer_detail::checkedAdd(
      static_cast<std::size_t>(sharedWignerBytes),readerPersistentBytes,"shared grid and primary reader bytes overflow");
  if(sharedResidentBytes>maxWorkingBytes-spoolValidationBytes)
    throw std::runtime_error("configured working-memory cap cannot hold the shared grid and primary spool reader");
  const std::uint64_t remainingAfterShared=maxWorkingBytes-spoolValidationBytes-sharedResidentBytes;
  if(basePerWorkerWithReserve>remainingAfterShared)
    throw std::runtime_error("configured working-memory cap cannot hold the shared grid and one base worker workspace");
  const std::size_t channelBatchSize=lmax<=32?4:
      (batchTwoScratch<=remainingAfterShared-basePerWorkerWithReserve?2:1);
  const std::size_t batchScratch=channelBatchSize>=4?batchFourScratch:(channelBatchSize==2?batchTwoScratch:0);
  if(batchScratch>std::numeric_limits<std::size_t>::max()-basePerWorkerBytes)
    throw std::overflow_error("same-spin batch scratch overflows addressable memory");
  const std::size_t knownPerWorkerBytes=basePerWorkerBytes+batchScratch;
  if(transformScratchBytes>std::numeric_limits<std::size_t>::max()-knownPerWorkerBytes)
    throw std::overflow_error("global transform scratch reserve overflows addressable memory");
  const std::size_t perWorkerBytes=knownPerWorkerBytes+static_cast<std::size_t>(transformScratchBytes);
  const std::size_t workersByMemory=perWorkerBytes==0?0:
      static_cast<std::size_t>((maxWorkingBytes-spoolValidationBytes-sharedResidentBytes)/perWorkerBytes);
  const std::size_t workers=std::min({options.workers,selectedSides.size(),static_cast<std::size_t>(48),workersByMemory});
  if(workers==0)throw std::runtime_error("configured working-memory cap cannot hold one complete-band side/frequency workspace");
  if(workers>48/options.gshThreads)
    throw std::runtime_error("effective global workers times inner GSH threads exceeds the 48-thread compute cap");
  const std::size_t innerThreadScratchBytesTotal=consumer_detail::checkedMultiply(innerThreadScratchBytesPerWorker,
      workers,"total inner-thread TLS bytes overflow");
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
    std::uint64_t transformsExecuted=0,zeroMixingSkipped=0;
  };
  std::vector<Worker> contexts(workers);
  Eigen::setNbThreads(1);
  const auto setupStarted=Clock::now();
  auto grid=std::make_shared<CompleteBandGrid>(lmax,options.storedWignerValues?
      GSHTrans::WignerValues::Stored():GSHTrans::WignerValues::Generated());
  for(auto &ctx:contexts) {
    ctx.reader=std::make_unique<PoleResponseSpoolReader>(spoolPath);
    ctx.plan=std::make_unique<GlobalGridPlan>(lmax,grid,channelBatchSize,options.gshThreads);
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
      if(records.audit.adjointMixing==Complex{}){++ctx.zeroMixingSkipped;continue;}
      ++ctx.transformsExecuted;
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
  std::uint64_t transformsExecuted=0,zeroMixingSkipped=0;
  std::uint64_t transformChannelsExecuted=0,fluidZeroTransformChannelsSkipped=0;
  for(const auto &ctx:contexts){
    workerElapsedSum+=ctx.sideSeconds;readSum+=ctx.recordReadSeconds;
    rotationSum+=ctx.plan->timings().coefficientRotationLayoutSeconds;
    inverseSum+=ctx.plan->timings().inverseGshFftSeconds;
    contractionSum+=ctx.plan->timings().pointKernelSeconds;
    accumulationSum+=ctx.frequencyAccumulationSeconds;writeSum+=ctx.outputWriteSeconds;
    transformsExecuted+=ctx.transformsExecuted;zeroMixingSkipped+=ctx.zeroMixingSkipped;
    transformChannelsExecuted+=ctx.plan->timings().transformChannelsExecuted;
    fluidZeroTransformChannelsSkipped+=ctx.plan->timings().fluidZeroTransformChannelsSkipped;
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
    <<",\"angular_transform_backend\":\"GSHTrans::TransformKernel::Loop with FFT longitude synthesis; Eigen GEMM is not used for angular synthesis\""
    <<",\"fftw_planner_policy\":\"FFTWpp::Estimate\""
    <<",\"wigner_storage_policy\":\""<<(options.storedWignerValues?"stored":"generated")<<"\",\"shared_wigner_table_bytes\":"<<sharedWignerBytes
    <<",\"stored_wigner_table_formula\":\"(2L+1)*(5*(L+1)^2-10)*sizeof(double), for n=-2..2 and All orders\""
    <<",\"serial_grid_and_context_setup_wall_s\":"<<setupSeconds
    <<",\"fftw_plan_note\":\"GSHTrans/FFTW transform plans are thread-local and initialize lazily on each worker's first inverse transform; that setup is included in that worker's first GSH inverse/FFT duration, not the serial grid-construction interval.\""
    <<",\"same_spin_channel_batch_size\":"<<channelBatchSize
    <<",\"gsh_threads_requested\":"<<options.gshThreads<<",\"gsh_execution_policy\":\""
    <<(options.gshThreads==1?"sequential":"GSHTrans::Execution::Parallel(explicit)")<<"\",\"openmp_available\":"
    <<(GSHTrans::OpenMP::Available?"true":"false")<<",\"openmp_dynamic_teams_disabled\":"
    <<((options.gshThreads>1&&GSHTrans::OpenMP::Available)?"true":"false")<<",\"configured_worker_count\":"<<options.workers
    <<",\"requested_thread_product\":"<<(options.workers*options.gshThreads)<<",\"effective_worker_count\":"<<workers
    <<",\"effective_thread_product_bound\":"<<(workers*options.gshThreads)<<",\"compute_thread_product_cap\":48"
    <<",\"generated_wigner_tls_bytes_per_inner_thread\":"<<generatedWignerTlsBytes
    <<",\"fft_input_output_tls_bytes_per_inner_thread_bound\":"<<fftBuffersTlsBytes
    <<",\"thread_local_fft_batch_slots_bound\":"<<threadLocalBatchStorage
    <<",\"known_tls_bytes_per_inner_thread_bound\":"<<knownTlsBytesPerInnerThread
    <<",\"inner_thread_tls_bytes_per_worker_bound\":"<<innerThreadScratchBytesPerWorker
    <<",\"inner_thread_tls_bytes_total_bound\":"<<innerThreadScratchBytesTotal
    <<",\"output_write_chunk_points\":"<<options.longitudeTile<<",\"workers\":"<<workers
    <<",\"max_radial_sides_cap\":"<<options.maxRadialSides<<",\"configured_working_memory_gib\":"<<options.maxWorkingGiB
    <<",\"known_working_bytes_per_worker\":"<<knownPerWorkerBytes
    <<",\"reader_persistent_metadata_bytes_per_instance\":"<<readerPersistentBytes
    <<",\"primary_and_worker_reader_count\":"<<(workers+1)
    <<",\"decoded_two_pole_record_bytes_per_worker\":"<<decodedRecordsBytes
    <<",\"two_frame_phase_norm_cache_bytes_per_worker\":"<<framePhaseCacheBytes
    <<",\"configured_transform_scratch_reserve_gib_per_worker\":"<<options.transformScratchGiBPerWorker
    <<",\"estimated_working_bytes_per_worker\":"<<perWorkerBytes
    <<",\"estimated_total_working_bytes\":"<<(sharedResidentBytes+static_cast<std::uint64_t>(workers)*perWorkerBytes+spoolValidationBytes)
    <<",\"spool_offset_validation_scratch_bytes\":"<<spoolValidationBytes
    <<",\"output_bytes\":"<<bytes<<",\"output_layout\":\"radial-side-major; theta-major angular prefix; point-major channels A,C,F,L,N,kappa; complete matching-frequency products accumulated per native side before output\""
    <<",\"kernel_is_full_earth_volume_integral\":false,\"source_contacts_included\":false"
    <<",\"binary_file\":\""<<std::filesystem::path(outputPath).filename().string()<<"\",\"elapsed_s\":"<<elapsed
    <<",\"timings\":{\"serial_grid_and_context_setup_wall_s\":"<<setupSeconds<<",\"worker_execution_wall_s\":"<<workerWallSeconds
    <<",\"frequency_side_transforms_executed\":"<<transformsExecuted<<",\"frequency_side_exact_zero_mixing_skipped\":"<<zeroMixingSkipped
    <<",\"transform_channels_executed\":"<<transformChannelsExecuted<<",\"fluid_exact_zero_channels_skipped\":"<<fluidZeroTransformChannelsSkipped
    <<",\"total_wall_s\":"<<elapsed<<",\"worker_elapsed_sum_s\":"<<workerElapsedSum
    <<",\"worker_stage_elapsed_sums\":{\"spool_reads\":"<<readSum<<",\"coefficient_rotation_layout\":"<<rotationSum
    <<",\"gsh_inverse_fft\":"<<inverseSum<<",\"point_kernel_contractions\":"<<contractionSum
    <<",\"frequency_accumulation\":"<<accumulationSum<<",\"output_writes\":"<<writeSum<<"},\"workers\":[";
  for(std::size_t i=0;i<contexts.size();++i){if(i)manifest<<',';const auto &ctx=contexts[i];const auto &t=ctx.plan->timings();
    manifest<<"{\"side_elapsed_s\":"<<ctx.sideSeconds<<",\"spool_reads_s\":"<<ctx.recordReadSeconds
      <<",\"coefficient_rotation_layout_s\":"<<t.coefficientRotationLayoutSeconds
      <<",\"gsh_inverse_fft_s\":"<<t.inverseGshFftSeconds<<",\"point_kernel_contractions_s\":"<<t.pointKernelSeconds
      <<",\"transform_channels_executed\":"<<t.transformChannelsExecuted<<",\"fluid_exact_zero_channels_skipped\":"<<t.fluidZeroTransformChannelsSkipped
      <<",\"frequency_accumulation_s\":"<<ctx.frequencyAccumulationSeconds<<",\"output_writes_s\":"<<ctx.outputWriteSeconds<<'}';}
  manifest<<"]},\"worker_timing_note\":\"Per-worker durations overlap in wall time and their sums are not a non-overlapping CPU-time measure. Total wall includes serial context setup and worker execution.\""
    <<",\"worker_owned_side_seconds\":[";
  for(std::size_t i=0;i<contexts.size();++i){if(i)manifest<<',';manifest<<contexts[i].sideSeconds;}
  manifest<<"],\"known_side_accumulator_bytes\":"<<fullSideAccumulatorBytes
    <<",\"memory_note\":\"Preflight charges a Stored Wigner table once as shared immutable grid memory, or Generated-mode angular scratch per worker; per-worker bounds include twelve complex fields, a frequency-result buffer, full-side accumulator, twelve coefficient arrays, two frame Wigners and phase/norm caches, chosen same-spin batch scratch, reader metadata, per-inner-thread Generated Wigner TLS and conservative FFT input/output buffers, and configured transform scratch reserve. FFTW/GSHTrans/Eigen/SEM/filesystem/thread-stack peak use must still be measured; this is not a peak-RSS guarantee.\"}\n";
  if(!manifest)throw std::runtime_error("failed writing global kernel manifest");
  ConsumerResult result;
  result.angularPointCount=sampleAngles;result.fullAngularPointCount=fullAngles;result.frequencyCount=metadata.frequencyCount();
  result.radialSideCount=selectedSides.size();result.outputBytes=bytes;result.elapsedSeconds=elapsed;
  result.serialGridContextSetupSeconds=setupSeconds;result.workerExecutionWallSeconds=workerWallSeconds;
  result.sharedWignerBytes=sharedWignerBytes;
  result.estimatedWorkingBytes=sharedResidentBytes+static_cast<std::uint64_t>(workers)*perWorkerBytes+spoolValidationBytes;
  result.frequencySideTransformsExecuted=transformsExecuted;result.frequencySideZeroMixingSkipped=zeroMixingSkipped;
  result.sameSpinBatchSize=channelBatchSize;result.gshThreads=options.gshThreads;result.workerCount=workers;
  result.requestedThreadProduct=options.workers*options.gshThreads;result.effectiveThreadProductBound=workers*options.gshThreads;
  result.innerThreadScratchBytesPerWorker=innerThreadScratchBytesPerWorker;
  result.innerThreadScratchBytesTotal=innerThreadScratchBytesTotal;
  return result;
}
} // namespace dspec_global_gsh
