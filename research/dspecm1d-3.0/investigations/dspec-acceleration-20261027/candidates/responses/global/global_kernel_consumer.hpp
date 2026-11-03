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

namespace dspec_global_gsh {
struct ConsumerOptions {
  std::size_t sampleAngleCount=0; // zero means the complete bandwidth grid
  std::size_t longitudeTile=128;
  std::size_t workers=1;
  double maxOutputGiB=8.0;
};
struct ConsumerResult {
  std::size_t angularPointCount=0,fullAngularPointCount=0,frequencyCount=0,radialSideCount=0;
  std::uint64_t outputBytes=0;
  double elapsedSeconds=0.0;
};

/** Replays compact source/unit-adjoint pole responses on the complete-band
 * angular grid (or a theta-major prefix such as the established first-64
 * diagnostic). Every output tile owns its full frequency sum and is written
 * once after the final matching-frequency product. */
inline ConsumerResult runPoleResponseSpoolKernel(
    const std::string &spoolPath,const Full1D::SEM &sem,
    const std::string &outputPath,ConsumerOptions options={}) {
  if(options.workers==0||options.workers>48||options.longitudeTile==0||
     !(std::isfinite(options.maxOutputGiB)&&options.maxOutputGiB>0.0))
    throw std::invalid_argument("invalid global kernel consumer worker, tile, or output limits");
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
  if(sides.empty()||sampleAngles>std::numeric_limits<std::size_t>::max()/sides.size()/6/sizeof(double))
    throw std::overflow_error("global kernel output size overflows size_t");
  const std::uint64_t bytes=static_cast<std::uint64_t>(sampleAngles)*sides.size()*6*sizeof(double);
  constexpr double gib=1024.0*1024.0*1024.0;
  if(options.maxOutputGiB>static_cast<double>(std::numeric_limits<std::uint64_t>::max())/gib||
     bytes>static_cast<std::uint64_t>(options.maxOutputGiB*gib))
    throw std::runtime_error("global kernel output exceeds configured byte cap");
  const auto parent=std::filesystem::path(outputPath).parent_path();
  if(!parent.empty())std::filesystem::create_directories(parent);
  if(std::filesystem::exists(outputPath)||std::filesystem::exists(outputPath+".json"))
    throw std::runtime_error("global kernel output or manifest path already exists");
  const auto space=std::filesystem::space(parent.empty()?std::filesystem::current_path():parent);
  if(bytes>space.available)throw std::runtime_error("insufficient free disk space for global kernel output");

  auto quadrature=GaussQuad::Quadrature1D<double>(GaussQuad::LegendrePolynomial<double>{}.GaussQuadrature(static_cast<int>(nTheta)));
  quadrature.Transform([](double x){return std::acos(-x);},[](double){return 1.0;});
  std::size_t tilesPerSide=0;
  for(std::size_t row=0;row<nTheta&&row*nPhi<sampleAngles;++row) {
    const std::size_t inRow=std::min(nPhi,sampleAngles-row*nPhi);
    tilesPerSide+=(inRow+options.longitudeTile-1)/options.longitudeTile;
  }
  if(tilesPerSide==0||sides.size()>std::numeric_limits<std::size_t>::max()/tilesPerSide)
    throw std::overflow_error("global kernel tile count overflows size_t");
  const std::size_t tileCount=sides.size()*tilesPerSide;
  const std::size_t workers=std::min(options.workers,tileCount);
  {
    std::ofstream initial(outputPath,std::ios::binary|std::ios::trunc);
    if(!initial)throw std::runtime_error("cannot create global kernel output");
    initial.seekp(static_cast<std::streamoff>(bytes-1));initial.put('\0');initial.close();
    if(!initial)throw std::runtime_error("cannot size global kernel output");
  }

  struct Worker {
    std::unique_ptr<PoleResponseSpoolReader> reader;
    std::unique_ptr<TilePlan> plan;
    std::fstream output;
    double seconds=0.0;
  };
  std::vector<Worker> contexts(workers);
  const auto started=std::chrono::steady_clock::now();
  runOwnedTiles(tileCount,workers,[&](std::size_t tile,std::size_t worker) {
    auto &ctx=contexts[worker];
    if(!ctx.reader) {
      ctx.reader=std::make_unique<PoleResponseSpoolReader>(spoolPath);
      TileLimits limits;limits.longitudeCount=options.longitudeTile;limits.maximumWorkers=workers;
      ctx.plan=std::make_unique<TilePlan>(lmax,limits);
      ctx.output.open(outputPath,std::ios::binary|std::ios::in|std::ios::out);
      if(!ctx.output)throw std::runtime_error("worker cannot open global kernel output");
    }
    const auto workStart=std::chrono::steady_clock::now();
    const std::size_t side=tile/tilesPerSide;
    std::size_t local=tile%tilesPerSide,row=0,phiBegin=0,pointCount=0;
    for(;row<nTheta;++row) {
      if(row*nPhi>=sampleAngles)throw std::logic_error("global tile maps beyond angular prefix");
      const std::size_t inRow=std::min(nPhi,sampleAngles-row*nPhi);
      const std::size_t rowTiles=(inRow+options.longitudeTile-1)/options.longitudeTile;
      if(local<rowTiles){phiBegin=local*options.longitudeTile;pointCount=std::min(options.longitudeTile,inRow-phiBegin);break;}
      local-=rowTiles;
    }
    if(row>=nTheta||pointCount==0)throw std::logic_error("global tile index is invalid");
    FrequencyTileAccumulator accumulator(pointCount);
    double maxBasisResidual=0.0;
    for(std::size_t bin=0;bin<ctx.reader->frequencyCount();++bin) {
      auto records=ctx.reader->read(bin,side);
      maxBasisResidual=std::max(maxBasisResidual,records.audit.basisReconstructionRelative);
      if(records.audit.basisReconstructionRelative>1e-12)
        throw std::runtime_error("receiver unit-basis reconstruction exceeds the 1e-12 gate");
      if(records.audit.adjointMixing==Complex{})continue;
      const auto frequencyKernel=ctx.plan->evaluate(records.forward,records.unitAdjoint,
          records.l0,records.audit.adjointMixing,sem,quadrature.Points()[row],phiBegin,nPhi,pointCount);
      accumulator.add(frequencyKernel);
    }
    const std::size_t angularBegin=row*nPhi+phiBegin;
    const std::uint64_t offset=(static_cast<std::uint64_t>(side)*sampleAngles+angularBegin)*6*sizeof(double);
    ctx.output.clear();ctx.output.seekp(static_cast<std::streamoff>(offset));
    const auto &values=accumulator.values();
    ctx.output.write(reinterpret_cast<const char*>(values.data()),static_cast<std::streamsize>(values.size()*sizeof(values.front())));
    if(!ctx.output)throw std::runtime_error("global kernel tile write failed");
    ctx.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-workStart).count();
  });

  const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  for(auto &ctx:contexts)if(ctx.output.is_open())ctx.output.close();
  std::size_t maxCachedWignerBytes=0;
  for(const auto &ctx:contexts)if(ctx.plan)maxCachedWignerBytes=std::max(maxCachedWignerBytes,ctx.plan->cachedWignerBytes());
  std::ofstream manifest(outputPath+".json");manifest<<std::setprecision(17);
  manifest<<"{\"mode\":\"compact_response_global_gsh_frequency_tile_kernel\",\"degree\":"<<lmax
    <<",\"angular_rule\":\"ForBand(L,2,2.0) nodes and FastFFTSize(4L+1) longitudes; band-L global GSH synthesis\""
    <<",\"n_theta\":"<<nTheta<<",\"n_phi\":"<<nPhi<<",\"full_angular_points\":"<<fullAngles
    <<",\"sampled_angular_points\":"<<sampleAngles<<",\"angular_subset_is_complete\":"<<(sampleAngles==fullAngles?"true":"false")
    <<",\"radial_side_count\":"<<sides.size()<<",\"frequency_count\":"<<metadata.frequencyCount()
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
    <<",\"angular_tile_longitudes\":"<<options.longitudeTile<<",\"workers\":"<<workers
    <<",\"configured_wigner_cache_bytes_per_worker\":"<<512u*1024u*1024u
    <<",\"known_cached_wigner_bytes_max_per_worker\":"<<maxCachedWignerBytes
    <<",\"output_bytes\":"<<bytes<<",\"output_layout\":\"radial-side-major; theta-major angular prefix; point-major channels A,C,F,L,N,kappa; frequency sum complete per output tile\""
    <<",\"kernel_is_full_earth_volume_integral\":false,\"source_contacts_included\":false"
    <<",\"binary_file\":\""<<std::filesystem::path(outputPath).filename().string()<<"\",\"elapsed_s\":"<<elapsed
    <<",\"worker_owned_tile_seconds\":[";
  for(std::size_t i=0;i<contexts.size();++i){if(i)manifest<<',';manifest<<contexts[i].seconds;}
  manifest<<"],\"known_tile_accumulator_bytes_max\":"<<options.longitudeTile*6*sizeof(double)
    <<",\"memory_note\":\"Known tile accumulator and configured Wigner cache bounds are separate from SEM, filesystem, thread-stack and Eigen/transform internals.\"}\n";
  if(!manifest)throw std::runtime_error("failed writing global kernel manifest");
  return {sampleAngles,fullAngles,metadata.frequencyCount(),sides.size(),bytes,elapsed};
}
} // namespace dspec_global_gsh
