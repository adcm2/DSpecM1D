#pragma once

#include "../engine_headers/Windowed3DKernels.h"
#include <GSHTrans/GaussLegendreGrid.hpp>
#include <GSHTrans/Indexing.hpp>
#include <GSHTrans/Wigner.hpp>
#include <GaussQuad/All>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <thread>
#include <exception>
#include <mutex>
#include <utility>
#include <vector>

namespace dspec_global_gsh {
using Complex = std::complex<double>;
using PoleRecord = SPARSESPEC::PoleAngularCoefficientRecord3D;
using Wigner = GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,GSHTrans::Single>;

struct L0Record {
  Complex forwardRr{}, forwardHorizontalTrace{};
  Complex adjointRr{}, adjointHorizontalTrace{};
};

/** Complete-band, FFT-backed global grid used for each frequency/native side.
 * The theta rule and longitude count match the accepted bounded prototype. */
class CompleteBandGrid final : public GSHTrans::SphericalGrid<double,GSHTrans::All,GSHTrans::All> {
  using Base=GSHTrans::SphericalGrid<double,GSHTrans::All,GSHTrans::All>;
  static std::pair<std::vector<double>,std::vector<double>> nodes(int lmax) {
    auto q=GaussQuad::Quadrature1D<double>(GaussQuad::LegendrePolynomial<double>{}.GaussQuadrature(2*lmax+1));
    q.Transform([](double x){return std::acos(-x);},[](double){return 1.0;});
    return {q.Points(),q.Weights()};
  }
public:
  explicit CompleteBandGrid(int lmax):CompleteBandGrid(lmax,nodes(lmax)){}
private:
  CompleteBandGrid(int lmax,std::pair<std::vector<double>,std::vector<double>> q)
    :Base(lmax,2,std::move(q.first),std::move(q.second),FFTWpp::Measure,
          GSHTrans::Chunking::Automatic(),GSHTrans::WignerValues::Generated(),
          GSHTrans::TransformKernel::Loop(),GSHTrans::FastFFTSize(4*lmax+1)){}
};

/** Rotate local pole coefficients into global order, apply the FFT-backed
 * complete-band inverse GSH transform for each spin channel, then contract
 * matching-frequency fields on the native side. */
class GlobalGridPlan {
public:
  explicit GlobalGridPlan(int lmax):lmax_(lmax),grid_(lmax) {
    if(lmax_<1||lmax_>1536)throw std::invalid_argument("global GSH degree is outside the supported range");
    nTheta_=static_cast<std::size_t>(2*lmax_+1);
    nPhi_=static_cast<std::size_t>(GSHTrans::FastFFTSize(4*lmax_+1));
    if(nTheta_>std::numeric_limits<std::size_t>::max()/nPhi_)throw std::overflow_error("global grid size overflows size_t");
    fieldSize_=nTheta_*nPhi_;
  }
  std::size_t fieldSize()const noexcept{return fieldSize_;}

  std::vector<std::array<double,6>> evaluate(const PoleRecord &forward,const PoleRecord &unitAdjoint,
      const L0Record &l0,Complex adjointMixing,const Full1D::SEM &sem) {
    if(forward.adjoint||!unitAdjoint.adjoint||forward.element!=unitAdjoint.element||
       forward.radiusModel!=unitAdjoint.radiusModel||forward.element<0||!finite(adjointMixing))
      throw std::invalid_argument("global GSH records or adjoint mixing do not match");
    std::array<const PoleRecord*,2> records{&forward,&unitAdjoint};
    std::array<std::vector<Complex>,12> fields;
    constexpr std::array<int,6> spins{0,0,1,-1,2,-2};
    for(int r=0;r<2;++r) {
      const auto &rec=*records[static_cast<std::size_t>(r)];
      const auto &rotation=rotationTable(rec.frame);
      const bool fluid=sem.elementIsFluid(rec.element);
      for(int component=0;component<6;++component) {
        const int spin=spins[static_cast<std::size_t>(component)];
        GSHTrans::GSHIndices<GSHTrans::All> indices(lmax_,lmax_,spin);
        std::vector<Complex> coefficients(static_cast<std::size_t>(indices.Size()),Complex{});
        for(int degree=std::abs(spin);degree<=lmax_;++degree) {
          const double norm=std::sqrt((2.0*degree+1.0)/(4.0*std::numbers::pi));
          for(int m=-degree;m<=degree;++m) {
            Complex value{};
            for(int mp=-std::min(2,degree);mp<=std::min(2,degree);++mp)
              value+=rotation.values->operator[](mp)[degree][m]*coefficient(rec,fluid,mp,degree,component)/norm;
            value*=std::exp(Complex(0.0,-static_cast<double>(m)*rotation.longitude));
            coefficients[static_cast<std::size_t>(indices.Index(degree,m))]=value*(r==1?adjointMixing:Complex{1.0,0.0});
          }
        }
        const std::size_t slot=static_cast<std::size_t>(r*6+component);
        fields[slot].resize(fieldSize_);
        grid_.InverseTransformation(lmax_,spin,coefficients,
          GSHTrans::Batch::Contiguous(1,static_cast<std::ptrdiff_t>(coefficients.size())),
          fields[slot],GSHTrans::Batch::Contiguous(1,static_cast<std::ptrdiff_t>(fieldSize_)),
          GSHTrans::Execution::Sequential());
      }
    }
    for(std::size_t p=0;p<fieldSize_;++p) {
      fields[0][p]+=l0.forwardRr;fields[1][p]+=l0.forwardHorizontalTrace;
      fields[6][p]+=l0.adjointRr*adjointMixing;fields[7][p]+=l0.adjointHorizontalTrace*adjointMixing;
    }
    std::vector<std::array<double,6>> result(fieldSize_);
    for(std::size_t p=0;p<fieldSize_;++p) {
      SPARSESPEC::Field3DSample sample;sample.radialElement=forward.element;
      sample.fluid=sem.elementIsFluid(forward.element);
      sample.forwardRawStrain={fields[0][p],fields[1][p],fields[2][p],fields[3][p],fields[4][p],fields[5][p]};
      sample.adjointDualStrain={fields[6][p],fields[7][p],fields[8][p],fields[9][p],fields[10][p],fields[11][p]};
      const auto kernel=SPARSESPEC::Windowed3DKernelEngine::pointKernel(sem,sample);
      result[p]={kernel.A,kernel.C,kernel.F,kernel.L,kernel.N,kernel.kappa};
      for(double v:result[p])if(!std::isfinite(v))throw std::runtime_error("global GSH kernel contains a nonfinite value");
    }
    return result;
  }
private:
  struct RotationTable {double longitude=0.0;std::unique_ptr<Wigner> values;};
  static bool finite(Complex z){return std::isfinite(z.real())&&std::isfinite(z.imag());}
  static std::pair<double,double> frameAngles(const SPARSESPEC::windowed3d_detail::PoleFrame3D &frame) {
    const auto &radial=frame.axes[2];return {std::acos(std::clamp(radial[2],-1.0,1.0)),std::atan2(radial[1],radial[0])};
  }
  const RotationTable &rotationTable(const SPARSESPEC::windowed3d_detail::PoleFrame3D &frame) {
    const auto [theta,longitude]=frameAngles(frame);
    const std::array<std::uint64_t,2> key{std::bit_cast<std::uint64_t>(theta),std::bit_cast<std::uint64_t>(longitude)};
    auto it=rotations_.find(key);if(it!=rotations_.end())return it->second;
    if(rotations_.size()>=2)throw std::runtime_error("global GSH frame cache limit exceeded");
    RotationTable value;value.longitude=longitude;value.values=std::make_unique<Wigner>(lmax_,lmax_,2,theta);
    return rotations_.emplace(key,std::move(value)).first->second;
  }
  int lmax_;CompleteBandGrid grid_;std::size_t nTheta_=0,nPhi_=0,fieldSize_=0;
  std::map<std::array<std::uint64_t,2>,RotationTable> rotations_;
};

/** Deterministic disjoint side ownership with one worker per independent
 * native radial side. Eigen's internal pool is pinned to one thread so this
 * outer pool is the sole parallel layer. Profiling counters must be collected
 * by the caller after joining workers. */
template<class Function>
void runOwnedTiles(std::size_t tileCount,std::size_t workers,Function &&function) {
  if(tileCount==0||workers==0||workers>48)
    throw std::invalid_argument("tile worker count must be between one and 48");
  workers=std::min(workers,tileCount);
  Eigen::setNbThreads(1);
  std::vector<std::thread> threads;threads.reserve(workers);
  std::mutex errorMutex;std::exception_ptr error;
  for(std::size_t worker=0;worker<workers;++worker) {
    const std::size_t begin=tileCount*worker/workers;
    const std::size_t end=tileCount*(worker+1)/workers;
    threads.emplace_back([&,worker,begin,end] {
      try {for(std::size_t tile=begin;tile<end;++tile)function(tile,worker);}
      catch(...) {std::lock_guard<std::mutex> lock(errorMutex);if(!error)error=std::current_exception();}
    });
  }
  for(auto &thread:threads)thread.join();
  if(error)std::rethrow_exception(error);
}

} // namespace dspec_global_gsh
