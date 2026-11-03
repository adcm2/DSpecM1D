#include "solvers/batch_response_cache.hpp"
#include <atomic>
#include <barrier>
#include <thread>
#include <iostream>
#include <vector>
using namespace dssolver;
void require(bool x){if(!x)throw std::runtime_error("cache concurrency assertion failed");}
int main(){
 constexpr int workers=8,keys=64,repeats=20;
 BatchResponseCache cache({1024,64u*1024u*1024u,1e-12});
 CacheDense rhs=CacheDense::Identity(4,4),basis=rhs;
 std::atomic<int> errors=0;
 std::barrier start(workers);
 std::vector<std::thread> threads;
 for(int worker=0;worker<workers;++worker)threads.emplace_back([&,worker]{
  try{
   start.arrive_and_wait();
   for(int repeat=0;repeat<repeats;++repeat)for(int i=0;i<keys;++i){
    auto key=BatchResponseKey::make(&cache,1,i,{double(i),-1.0});
    auto out=cache.getOrCompute(key,rhs,basis,[](const auto &f,const auto &b){
     std::this_thread::yield();return std::pair<CacheDense,CacheDense>{f*2.0,b*3.0};});
    require((out.forward-rhs*2.0).norm()==0.0);
    auto weighted=cache.deriveWeightedAdjoint(key,basis,basis*5.0);
    require((weighted.weightedAdjoint-basis*15.0).norm()==0.0);
   }
  }catch(...){++errors;}
 });
 for(auto &thread:threads)thread.join();threads.clear();
 require(errors==0&&cache.entryCount()==keys&&cache.retainedBytes()<=64u*1024u*1024u);
 auto key=BatchResponseKey::make(&cache,1,0,{0,-1});
 bool mismatch=false;try{cache.getOrCompute(key,rhs*7.0,basis,[](auto &f,auto &b){return std::pair{f,b};});}catch(const BatchResponseCache::BasisMismatch&){mismatch=true;}require(mismatch);
 // A callback exception must leave the cache unlocked and without an entry.
 auto errorKey=BatchResponseKey::make(&cache,1,999,{2,-1});
 bool caught=false;try{cache.getOrCompute(errorKey,rhs,basis,[](const auto&,const auto&)->std::pair<CacheDense,CacheDense>{throw std::runtime_error("injected");});}catch(const std::runtime_error&){caught=true;}
 require(caught&&cache.entryCount()==keys);
 cache.getOrCompute(errorKey,rhs,basis,[](auto &f,auto &b){return std::pair{f,b};});
 // Concurrent publication must enforce the aggregate entry cap.
 BatchResponseCache capped({1,1u<<20,1e-12});std::atomic<int> accepted=0,rejected=0;
 std::barrier capStart(workers);
 for(int i=0;i<workers;++i)threads.emplace_back([&,i]{capStart.arrive_and_wait();try{
  auto key=BatchResponseKey::make(&capped,1,i,{double(i),-1});
  capped.getOrCompute(key,rhs,basis,[](auto &f,auto &b){std::this_thread::yield();return std::pair{f,b};});++accepted;
 }catch(const BatchResponseCache::CapacityExceeded&){++rejected;}catch(...){++errors;}});
 for(auto &thread:threads)thread.join();require(accepted==1&&rejected==workers-1&&errors==0&&capped.entryCount()==1);
 std::cout<<"{\"passed\":true,\"workers\":8,\"response_requests\":"<<workers*keys*repeats<<",\"duplicate_key_publication\":true,\"weighted_reads_during_insert\":true,\"signature_rejection\":true,\"callback_exception_recovery\":true,\"concurrent_capacity_enforced\":true}\n";
}
