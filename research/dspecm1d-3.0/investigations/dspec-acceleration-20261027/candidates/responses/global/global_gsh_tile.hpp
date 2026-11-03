#pragma once

#include "../engine_headers/Windowed3DKernels.h"
#include <GSHTrans/Wigner.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
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

struct TileLimits {
  std::size_t longitudeCount=128;
  std::size_t maximumCachedWignerBytes=512u*1024u*1024u;
  std::size_t maximumCachedFrames=2;
  std::size_t maximumWorkers=48;
};

/**
 * Bounded global-frame GSH synthesis for one native radial side and one
 * theta/longitude tile. Local pole records retain the established spin and
 * frame convention. Only one degree's angular basis and coefficient blocks
 * are live at a time; Eigen performs the degree-block products internally.
 */
class TilePlan {
public:
  explicit TilePlan(int lmax,TileLimits limits={}) : lmax_(lmax),limits_(limits) {
    if(lmax_<1||lmax_>1536||limits_.longitudeCount==0||
       limits_.maximumCachedFrames==0||limits_.maximumCachedWignerBytes==0||
       limits_.maximumWorkers==0||limits_.maximumWorkers>48)
      throw std::invalid_argument("invalid global GSH tile limits");
    thetaBytes_=SPARSESPEC::Windowed3DFieldEngine::angularWignerBytesPerAngle(lmax_,2);
  }

  std::size_t longitudeCount() const noexcept { return limits_.longitudeCount; }
  std::size_t cachedWignerBytes()const noexcept{return cachedBytes_;}

  std::vector<std::array<double,6>> evaluate(
      const PoleRecord &forward,const PoleRecord &unitAdjoint,
      const L0Record &l0,double adjointMixing,const Full1D::SEM &sem,
      double theta,std::size_t phiBegin,std::size_t phiTotal,std::size_t requestedCount=0) {
    if(forward.adjoint||!unitAdjoint.adjoint||forward.element!=unitAdjoint.element||
       forward.radiusModel!=unitAdjoint.radiusModel||forward.element<0||
       phiTotal==0||phiBegin>=phiTotal||!std::isfinite(adjointMixing))
      throw std::invalid_argument("global GSH records or tile coordinates do not match");
    const std::size_t count=std::min({limits_.longitudeCount,phiTotal-phiBegin,
        requestedCount==0?limits_.longitudeCount:requestedCount});
    const auto &globalWigner=thetaTable(theta);
    const auto &forwardRotation=rotationTable(forward.frame);
    const auto &adjointRotation=rotationTable(unitAdjoint.frame);
    const double longitudeStep=2.0*std::acos(-1.0)/static_cast<double>(phiTotal);
    std::vector<double> phi(count);
    for(std::size_t p=0;p<count;++p)phi[p]=longitudeStep*static_cast<double>(phiBegin+p);

    std::array<std::vector<Complex>,12> fields;
    for(auto &field:fields)field.assign(count,Complex{});
    constexpr std::array<int,6> spin{0,0,1,-1,2,-2};
    for(int degree=0;degree<=lmax_;++degree) {
      std::array<std::unique_ptr<Eigen::MatrixXcd>,5> bases;
      for(int component=0;component<6;++component) {
        if(degree<std::abs(spin[static_cast<std::size_t>(component)]))continue;
        const int n=2*degree+1;
        const int basisSlot=spin[static_cast<std::size_t>(component)]+2;
        if(!bases[static_cast<std::size_t>(basisSlot)]) {
          auto basis=std::make_unique<Eigen::MatrixXcd>(static_cast<Eigen::Index>(count),n);
          for(int m=-degree;m<=degree;++m) {
            const Complex phase=std::exp(Complex(0.0,static_cast<double>(m)*phi.front()));
            Complex p=phase;
            const Complex step=std::exp(Complex(0.0,static_cast<double>(m)*longitudeStep));
            const Complex angular=globalWigner[spin[static_cast<std::size_t>(component)]][degree][m];
            for(std::size_t j=0;j<count;++j) {
              (*basis)(static_cast<Eigen::Index>(j),m+degree)=angular*p;
              p*=step;
            }
          }
          bases[static_cast<std::size_t>(basisSlot)]=std::move(basis);
        }
        Eigen::MatrixXcd fcoeff(n,1),acoeff(n,1);
        const bool fluid=sem.elementIsFluid(forward.element);
        rotateDegree(forward,fluid,forwardRotation,degree,component,fcoeff);
        rotateDegree(unitAdjoint,fluid,adjointRotation,degree,component,acoeff);
        acoeff*=adjointMixing;
        const Eigen::MatrixXcd ffield=(*bases[static_cast<std::size_t>(basisSlot)])*fcoeff;
        const Eigen::MatrixXcd afield=(*bases[static_cast<std::size_t>(basisSlot)])*acoeff;
        const std::size_t fslot=static_cast<std::size_t>(component);
        const std::size_t aslot=6+fslot;
        for(std::size_t j=0;j<count;++j) {
          fields[fslot][j]+=ffield(static_cast<Eigen::Index>(j),0);
          fields[aslot][j]+=afield(static_cast<Eigen::Index>(j),0);
        }
      }
    }
    for(std::size_t j=0;j<count;++j) {
      fields[0][j]+=l0.forwardRr;
      fields[1][j]+=l0.forwardHorizontalTrace;
      fields[6][j]+=l0.adjointRr*adjointMixing;
      fields[7][j]+=l0.adjointHorizontalTrace*adjointMixing;
    }

    std::vector<std::array<double,6>> result(count);
    for(std::size_t j=0;j<count;++j) {
      SPARSESPEC::Field3DSample sample;
      sample.radialElement=forward.element;
      sample.fluid=sem.elementIsFluid(forward.element);
      sample.forwardRawStrain={fields[0][j],fields[1][j],fields[2][j],fields[3][j],fields[4][j],fields[5][j]};
      sample.adjointDualStrain={fields[6][j],fields[7][j],fields[8][j],fields[9][j],fields[10][j],fields[11][j]};
      const auto kernel=SPARSESPEC::Windowed3DKernelEngine::pointKernel(sem,sample);
      result[j]={kernel.A,kernel.C,kernel.F,kernel.L,kernel.N,kernel.kappa};
      for(double value:result[j])if(!std::isfinite(value))
        throw std::runtime_error("global GSH kernel tile contains a nonfinite value");
    }
    return result;
  }

private:
  using FrameKey=std::array<std::uint64_t,2>;
  struct RotationTable { double longitude=0.0;std::unique_ptr<Wigner> values;std::size_t bytes=0; };

  static std::pair<double,double> frameAngles(const SPARSESPEC::windowed3d_detail::PoleFrame3D &frame) {
    const auto &radial=frame.axes[2];
    return {std::acos(std::clamp(radial[2],-1.0,1.0)),std::atan2(radial[1],radial[0])};
  }

  const Wigner &thetaTable(double theta) {
    const auto key=std::bit_cast<std::uint64_t>(theta);
    if(thetaTable_&&thetaKey_==key)return *thetaTable_;
    if(thetaTable_){thetaTable_.reset();cachedBytes_-=thetaBytes_;}
    ensureCache(thetaBytes_);
    thetaTable_=std::make_unique<Wigner>(lmax_,lmax_,2,theta);
    thetaKey_=key;
    cachedBytes_+=thetaBytes_;
    return *thetaTable_;
  }

  const RotationTable &rotationTable(const SPARSESPEC::windowed3d_detail::PoleFrame3D &frame) {
    const auto [theta,longitude]=frameAngles(frame);
    const FrameKey key{std::bit_cast<std::uint64_t>(theta),std::bit_cast<std::uint64_t>(longitude)};
    auto it=rotationTables_.find(key);
    if(it!=rotationTables_.end())return it->second;
    if(rotationTables_.size()>=limits_.maximumCachedFrames)
      throw std::runtime_error("global GSH rotation-frame cache limit exceeded");
    ensureCache(thetaBytes_);
    RotationTable table;table.longitude=longitude;table.bytes=thetaBytes_;
    table.values=std::make_unique<Wigner>(lmax_,lmax_,2,theta);
    cachedBytes_+=table.bytes;
    return rotationTables_.emplace(key,std::move(table)).first->second;
  }

  void ensureCache(std::size_t extra) const {
    if(extra>limits_.maximumCachedWignerBytes-cachedBytes_)
      throw std::runtime_error("global GSH Wigner cache exceeds its configured byte limit");
  }

  static Complex coefficient(const PoleRecord &record,bool fluid,int localOrder,int degree,int component) {
    const auto it=std::lower_bound(record.labels.begin(),record.labels.end(),std::pair{localOrder,degree});
    if(it==record.labels.end()||*it!=std::pair{localOrder,degree})return {};
    const std::size_t index=static_cast<std::size_t>(it-record.labels.begin());
    const auto &v=record.rawStrain[index];
    if(fluid)return component==0?v[0]+v[1]:Complex{};
    return v[static_cast<std::size_t>(component)];
  }

  void rotateDegree(const PoleRecord &record,bool fluid,const RotationTable &rotation,
                    int degree,int component,Eigen::MatrixXcd &out) const {
    const double norm=std::sqrt((2.0*degree+1.0)/(4.0*std::acos(-1.0)));
    Eigen::MatrixXcd local=Eigen::MatrixXcd::Zero(5,1);
    for(int mp=-2;mp<=2;++mp)
      if(std::abs(mp)<=degree)local(mp+2,0)=coefficient(record,fluid,mp,degree,component);
    for(int m=-degree;m<=degree;++m) {
      Complex value{};
      for(int mp=-2;mp<=2;++mp)
        if(std::abs(mp)<=degree)value+=rotation.values->operator[](mp)[degree][m]*local(mp+2,0)/norm;
      value*=std::exp(Complex(0.0,-static_cast<double>(m)*rotation.longitude));
      out(m+degree,0)=value;
    }
  }

  int lmax_;
  TileLimits limits_;
  std::size_t thetaBytes_=0;
  mutable std::size_t cachedBytes_=0;
  std::uint64_t thetaKey_=0;
  std::unique_ptr<Wigner> thetaTable_;
  std::map<FrameKey,RotationTable> rotationTables_;
};

/** Frequency summation owned by one output tile; callers write only after the
 * final matching-frequency product has been accumulated. */
class FrequencyTileAccumulator {
public:
  explicit FrequencyTileAccumulator(std::size_t points) : sums_(points) {
    if(points==0)throw std::invalid_argument("frequency tile must contain points");
  }
  void add(const std::vector<std::array<double,6>> &frequencyKernel) {
    if(frequencyKernel.size()!=sums_.size())throw std::invalid_argument("frequency tile size changed");
    for(std::size_t p=0;p<sums_.size();++p)
      for(std::size_t c=0;c<6;++c)sums_[p][c]+=frequencyKernel[p][c];
  }
  const std::vector<std::array<double,6>> &values() const noexcept { return sums_; }
private:
  std::vector<std::array<double,6>> sums_;
};

/** Deterministic disjoint tile ownership with one worker per independent
 * angular tile range. Eigen's internal pool is pinned to one thread so this
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
