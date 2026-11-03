#pragma once

// Exact, bounded response reuse for one explicitly scoped batch request.
// This cache stores solved response blocks, never sparse matrices or factors.
#include <Eigen/Core>
#include <Eigen/QR>

#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dssolver {
using CacheComplex = std::complex<double>;
using CacheDense = Eigen::MatrixXcd;

struct BatchResponseKey {
  std::uintptr_t semIdentity=0;
  int family=0;
  int degree=0;
  std::uint64_t zRealBits=0;
  std::uint64_t zImagBits=0;

  static BatchResponseKey make(const void *sem,int family,int degree,CacheComplex z) {
    if(!sem||!std::isfinite(z.real())||!std::isfinite(z.imag()))
      throw std::invalid_argument("batch response key requires a SEM identity and finite frequency");
    return {reinterpret_cast<std::uintptr_t>(sem),family,degree,
            std::bit_cast<std::uint64_t>(z.real()),std::bit_cast<std::uint64_t>(z.imag())};
  }
  friend bool operator<(const BatchResponseKey &a,const BatchResponseKey &b) {
    if(a.semIdentity!=b.semIdentity)return a.semIdentity<b.semIdentity;
    if(a.family!=b.family)return a.family<b.family;
    if(a.degree!=b.degree)return a.degree<b.degree;
    if(a.zRealBits!=b.zRealBits)return a.zRealBits<b.zRealBits;
    return a.zImagBits<b.zImagBits;
  }
};

/**
 * Scoped memoization for repeated-objective adjoint batches. On a cache miss,
 * `computeResponses(forwardRhs, unitAdjointRhs)` must use one factorization to
 * return the forward response and Hermitian-adjoint response to every column
 * of the unit basis. On a hit, the global weighted adjoint RHS is represented
 * in that basis and its response is formed by the same thin mixing matrix.
 *
 * A key collision with different forward or unit-basis RHS contents, a rank
 * deficient basis, or a failed reconstruction is a hard error. Callers may
 * catch `BasisMismatch` and explicitly run the ordinary weighted adjoint.
 * The object is intended to be stack-scoped to one batch request; destruction
 * releases every response and signature. Entry and byte caps bound retained
 * state even if a caller exceeds the documented request limits.
 */
class BatchResponseCache {
public:
  struct Result {
    CacheDense forward;
    CacheDense unitAdjoint;
    bool cacheHit=false;
  };
  struct WeightedResult {
    CacheDense mixing;
    CacheDense weightedAdjoint;
    double basisReconstructionRelative=0.0;
  };
  struct Limits {
    std::size_t maximumEntries=2048;
    std::size_t maximumBytes=512u*1024u*1024u;
    double maximumBasisReconstructionRelative=1e-12;
  };
  struct BasisMismatch : std::runtime_error { using std::runtime_error::runtime_error; };
  struct CapacityExceeded : std::runtime_error { using std::runtime_error::runtime_error; };

  BatchResponseCache() : BatchResponseCache(Limits{}) {}
  explicit BatchResponseCache(Limits limits) : limits_(limits) {
    if(limits_.maximumEntries==0||limits_.maximumBytes==0||
       !(limits_.maximumBasisReconstructionRelative>0.0)||
       !std::isfinite(limits_.maximumBasisReconstructionRelative))
      throw std::invalid_argument("invalid batch response cache limits");
  }
  BatchResponseCache(const BatchResponseCache&)=delete;
  BatchResponseCache& operator=(const BatchResponseCache&)=delete;
  BatchResponseCache(BatchResponseCache&&)=default;
  BatchResponseCache& operator=(BatchResponseCache&&)=default;

  template<class ComputeResponses>
  Result getOrCompute(const BatchResponseKey &key,
                      const CacheDense &forwardRhs,
                      const CacheDense &unitAdjointRhs,
                      ComputeResponses &&computeResponses) {
    validateRhs(forwardRhs,"forward RHS");
    validateRhs(unitAdjointRhs,"unit-basis adjoint RHS");
    if(unitAdjointRhs.rows()!=forwardRhs.rows())
      throw std::invalid_argument("batch response RHS row dimensions do not match");

    auto found=entries_.find(key);
    if(found!=entries_.end()) {
      const Entry &entry=found->second;
      if(!sameBits(entry.forwardRhsSignature,forwardRhs)||
         !sameBits(entry.unitAdjointRhsSignature,unitAdjointRhs))
        throw BasisMismatch("batch cache key reused with a different forward or unit-basis RHS");
      ++cacheHits_;
      Result result;
      result.forward=entry.forwardResponse;
      result.unitAdjoint=entry.unitAdjointResponse;
      result.cacheHit=true;
      return result;
    }

    ++cacheMisses_;
    if(entries_.size()>=limits_.maximumEntries)
      throw CapacityExceeded("batch response cache entry cap reached");
    auto responses=std::invoke(std::forward<ComputeResponses>(computeResponses),
                               forwardRhs,unitAdjointRhs);
    CacheDense forwardResponse=std::move(responses.first);
    CacheDense unitAdjointResponse=std::move(responses.second);
    if(forwardResponse.rows()!=forwardRhs.rows()||
       forwardResponse.cols()!=forwardRhs.cols()||
       unitAdjointResponse.rows()!=unitAdjointRhs.rows()||
       unitAdjointResponse.cols()!=unitAdjointRhs.cols()||
       !forwardResponse.allFinite()||!unitAdjointResponse.allFinite())
      throw std::runtime_error("computed batch response dimensions or values are invalid");

    Entry entry;
    entry.forwardRhsSignature=bitSignature(forwardRhs);
    entry.unitAdjointRhsSignature=bitSignature(unitAdjointRhs);
    entry.forwardResponse=std::move(forwardResponse);
    entry.unitAdjointResponse=std::move(unitAdjointResponse);
    entry.bytes=matrixBytes(forwardRhs)+matrixBytes(unitAdjointRhs)+
                matrixBytes(entry.forwardResponse)+matrixBytes(entry.unitAdjointResponse)+
                (entry.forwardRhsSignature.size()+entry.unitAdjointRhsSignature.size())*sizeof(std::uint64_t);
    if(entry.bytes>limits_.maximumBytes-currentBytes_)
      throw CapacityExceeded("batch response cache byte cap reached");
    currentBytes_+=entry.bytes;
    entries_.emplace(key,std::move(entry));

    Result result;
    result.forward=entries_.at(key).forwardResponse;
    result.unitAdjoint=entries_.at(key).unitAdjointResponse;
    result.cacheHit=false;
    return result;
  }

  WeightedResult deriveWeightedAdjoint(const BatchResponseKey &key,
                                       const CacheDense &unitAdjointRhs,
                                       const CacheDense &weightedAdjointRhs) const {
    validateRhs(unitAdjointRhs,"unit-basis adjoint RHS");
    validateRhs(weightedAdjointRhs,"weighted adjoint RHS");
    const auto found=entries_.find(key);
    if(found==entries_.end())throw std::out_of_range("no cached response for batch key");
    const Entry &entry=found->second;
    if(!sameBits(entry.unitAdjointRhsSignature,unitAdjointRhs))
      throw BasisMismatch("weighted adjoint uses a different unit-basis RHS signature");
    if(weightedAdjointRhs.rows()!=unitAdjointRhs.rows())
      throw std::invalid_argument("weighted adjoint RHS row dimension does not match its basis");
    WeightedResult result;
    result.mixing=representInBasis(unitAdjointRhs,weightedAdjointRhs);
    result.basisReconstructionRelative=
        reconstructionRelative(unitAdjointRhs,result.mixing,weightedAdjointRhs);
    if(!(result.basisReconstructionRelative<=limits_.maximumBasisReconstructionRelative))
      throw BasisMismatch("weighted adjoint RHS is outside the unit-basis span");
    maximumBasisReconstructionRelative_=std::max(
        maximumBasisReconstructionRelative_,result.basisReconstructionRelative);
    result.weightedAdjoint=entry.unitAdjointResponse*result.mixing;
    if(!result.weightedAdjoint.allFinite())
      throw std::runtime_error("derived weighted adjoint response is nonfinite");
    return result;
  }

  std::size_t entryCount() const noexcept { return entries_.size(); }
  std::size_t retainedBytes() const noexcept { return currentBytes_; }
  std::size_t cacheHits() const noexcept { return cacheHits_; }
  std::size_t cacheMisses() const noexcept { return cacheMisses_; }
  double maximumBasisReconstructionRelative() const noexcept {
    return maximumBasisReconstructionRelative_;
  }
  void clear() noexcept { entries_.clear();currentBytes_=0; }

private:
  struct Entry {
    std::vector<std::uint64_t> forwardRhsSignature;
    std::vector<std::uint64_t> unitAdjointRhsSignature;
    CacheDense forwardResponse;
    CacheDense unitAdjointResponse;
    std::size_t bytes=0;
  };
  static void validateRhs(const CacheDense &rhs,const char *name) {
    if(rhs.rows()<=0||rhs.cols()<=0||!rhs.allFinite())
      throw std::invalid_argument(std::string("invalid ")+name);
  }
  static std::size_t matrixBytes(const CacheDense &matrix) {
    const auto count=static_cast<std::size_t>(matrix.size());
    if(count>std::numeric_limits<std::size_t>::max()/sizeof(CacheComplex))
      throw CapacityExceeded("batch response matrix byte count overflows size_t");
    return count*sizeof(CacheComplex);
  }
  static std::vector<std::uint64_t> bitSignature(const CacheDense &matrix) {
    std::vector<std::uint64_t> bits;
    bits.reserve(static_cast<std::size_t>(2*matrix.size()+2));
    bits.push_back(static_cast<std::uint64_t>(matrix.rows()));
    bits.push_back(static_cast<std::uint64_t>(matrix.cols()));
    for(Eigen::Index j=0;j<matrix.cols();++j)
      for(Eigen::Index i=0;i<matrix.rows();++i) {
        bits.push_back(std::bit_cast<std::uint64_t>(matrix(i,j).real()));
        bits.push_back(std::bit_cast<std::uint64_t>(matrix(i,j).imag()));
      }
    return bits;
  }
  static bool sameBits(const std::vector<std::uint64_t> &signature,
                       const CacheDense &matrix) {
    const auto expected=bitSignature(matrix);
    return signature==expected;
  }
  static CacheDense representInBasis(const CacheDense &basis,const CacheDense &rhs) {
    Eigen::ColPivHouseholderQR<CacheDense> qr(basis);
    if(qr.rank()!=basis.cols())
      throw BasisMismatch("unit adjoint RHS basis is rank deficient");
    CacheDense mix=qr.solve(rhs);
    if(!mix.allFinite())throw BasisMismatch("unit adjoint basis solve is nonfinite");
    return mix;
  }
  static double reconstructionRelative(const CacheDense &basis,const CacheDense &mix,
                                       const CacheDense &rhs) {
    const double numerator=(basis*mix-rhs).norm();
    const double denominator=rhs.norm();
    if(denominator==0.0)
      return numerator==0.0?0.0:std::numeric_limits<double>::infinity();
    return numerator/denominator;
  }

  Limits limits_;
  std::map<BatchResponseKey,Entry> entries_;
  std::size_t currentBytes_=0;
  std::size_t cacheHits_=0,cacheMisses_=0;
  mutable double maximumBasisReconstructionRelative_=0.0;
};

} // namespace dssolver
