#include "batch_response_cache.hpp"
#include <cassert>
#include <cmath>
#include <stdexcept>

int main() {
  int semToken=0;
  const auto key=dssolver::BatchResponseKey::make(&semToken,1,2,{3.0,-0.125});
  dssolver::BatchResponseCache cache;
  dssolver::CacheDense f=dssolver::CacheDense::Zero(4,1);f(0,0)=1.0;
  dssolver::CacheDense basis=dssolver::CacheDense::Zero(4,2);
  basis(0,0)=1.0;basis(1,1)=1.0;
  dssolver::CacheDense mixing(2,2);
  mixing(0,0)={0.3,0.1};mixing(0,1)={-0.2,0.4};
  mixing(1,0)={0.5,-0.2};mixing(1,1)={0.7,0.0};
  const dssolver::CacheDense weighted=basis*mixing;
  dssolver::CacheDense forwardResponse(4,1);forwardResponse<<2.0,3.0,4.0,5.0;
  dssolver::CacheDense basisResponse(4,2);basisResponse<<1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0;
  const auto first=cache.getOrCompute(key,f,basis,[&](const auto &,const auto &) {
    return std::pair{forwardResponse,basisResponse};
  });
  assert(!first.cacheHit&&cache.entryCount()==1);
  const auto derived=cache.deriveWeightedAdjoint(key,basis,weighted);
  assert(derived.basisReconstructionRelative<1e-14);
  assert((derived.weightedAdjoint-basisResponse*mixing).norm()<1e-13);
  const auto second=cache.getOrCompute(key,f,basis,[&](const auto &,const auto &) -> std::pair<dssolver::CacheDense,dssolver::CacheDense> {
    throw std::runtime_error("cache miss unexpectedly invoked callback");
  });
  assert(second.cacheHit);
  assert(cache.cacheMisses()==1&&cache.cacheHits()==1);
  assert(cache.maximumBasisReconstructionRelative()<1e-14);
  auto changedForward=f;changedForward(1,0)=1.0;
  bool mismatch=false;
  try { (void)cache.getOrCompute(key,changedForward,basis,[&](const auto &,const auto &) {
    return std::pair{forwardResponse,basisResponse};
  }); } catch(const dssolver::BatchResponseCache::BasisMismatch &) { mismatch=true; }
  assert(mismatch);
  auto outOfSpan=weighted;outOfSpan(3,0)=1.0;
  mismatch=false;
  try { (void)cache.deriveWeightedAdjoint(key,basis,outOfSpan); }
  catch(const dssolver::BatchResponseCache::BasisMismatch &) { mismatch=true; }
  assert(mismatch);
  return 0;
}
