#pragma once
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
namespace mathprofile {
using Clock=std::chrono::steady_clock;
struct Entry {std::uint64_t inclusive=0,exclusive=0,calls=0;};
struct Frame {const char*name=nullptr;Clock::time_point start;std::uint64_t children=0;};
inline thread_local bool active=false;
inline thread_local std::array<Frame,128> stack;
inline thread_local std::size_t depth=0;
inline thread_local std::map<std::string,Entry> entries;
inline thread_local std::map<std::string,std::uint64_t> counters;
// Worker stage times are sums of elapsed worker durations, not wall times.
// Save/restore also handles the OpenMP master, which already has open scopes.
struct WorkerData {
 std::map<std::string,Entry> entries;
 std::map<std::string,std::uint64_t> counters;
};
struct Snapshot {
 std::map<std::string,Entry> entries;
 std::map<std::string,std::uint64_t> counters;
};
struct WorkerCapture {
 bool previousActive;std::size_t previousDepth;
 std::array<Frame,128> previousStack;
 WorkerData previous;WorkerData &output;
 WorkerCapture(bool enabled,WorkerData &data):previousActive(active),previousDepth(depth),
     previousStack(stack),output(data){
  previous.entries.swap(entries);previous.counters.swap(counters);depth=0;active=enabled;
 }
 ~WorkerCapture(){
  if(depth)std::terminate();
  output.entries.swap(entries);output.counters.swap(counters);
  previous.entries.swap(entries);previous.counters.swap(counters);
  stack=previousStack;depth=previousDepth;active=previousActive;
 }
 WorkerCapture(const WorkerCapture&)=delete;
};
inline void mergeWorker(const WorkerData &data){
 for(const auto &[name,value]:data.entries){auto &e=entries["parallel_worker."+name];
  e.inclusive+=value.inclusive;e.exclusive+=value.exclusive;e.calls+=value.calls;}
 for(const auto &[name,value]:data.counters)counters[name]+=value;
}
#ifdef MATHPROFILE_DISABLED
struct Scope { explicit Scope(const char*) noexcept {} void stop() noexcept {} };
inline void count(const char*,std::uint64_t=1) noexcept {}
inline Snapshot snapshot(){return {};}
inline Snapshot difference(const Snapshot&,const Snapshot&){return {};}
inline void accumulate(Snapshot&,const Snapshot&) noexcept {}
struct Root { explicit Root(const char*,int=0) noexcept {} };
#else
struct Scope {
 bool enabled,finished=false;std::size_t index=0;
 explicit Scope(const char*n):enabled(active){if(enabled){index=depth++;if(depth>stack.size())std::terminate();stack[index]={n,Clock::now(),0};}}
 ~Scope(){stop();}
 void stop() noexcept {if(!enabled||finished)return;if(depth!=index+1)std::terminate();const auto end=Clock::now();auto &f=stack[index];const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(end-f.start).count();--depth;auto &e=entries[f.name];e.inclusive+=elapsed;e.exclusive+=elapsed-f.children;++e.calls;if(depth)stack[depth-1].children+=elapsed;finished=true;}
 Scope(const Scope&)=delete;
};
inline void count(const char*n,std::uint64_t v=1){if(active)counters[n]+=v;}
inline Snapshot snapshot(){return {entries,counters};}
inline Snapshot difference(const Snapshot &current,const Snapshot &previous){
 Snapshot out;
 for(const auto &[name,value]:current.entries){
  const auto found=previous.entries.find(name);
  const Entry old=found==previous.entries.end()?Entry{}:found->second;
  if(value.inclusive<old.inclusive||value.exclusive<old.exclusive||value.calls<old.calls)std::terminate();
  const Entry delta{value.inclusive-old.inclusive,value.exclusive-old.exclusive,value.calls-old.calls};
  if(delta.inclusive||delta.exclusive||delta.calls)out.entries.emplace(name,delta);
 }
 for(const auto &[name,value]:current.counters){
  const auto found=previous.counters.find(name);const std::uint64_t old=found==previous.counters.end()?0:found->second;
  if(value<old)std::terminate();const auto delta=value-old;if(delta)out.counters.emplace(name,delta);
 }
 for(const auto &[name,value]:previous.entries)if(!current.entries.contains(name))std::terminate();
 for(const auto &[name,value]:previous.counters)if(!current.counters.contains(name))std::terminate();
 return out;
}
inline void accumulate(Snapshot &target,const Snapshot &delta){
 for(const auto &[name,value]:delta.entries){auto &entry=target.entries[name];entry.inclusive+=value.inclusive;entry.exclusive+=value.exclusive;entry.calls+=value.calls;}
 for(const auto &[name,value]:delta.counters)target.counters[name]+=value;
}
struct Root {
 const char*solver;int objective;
 explicit Root(const char*s,int o=0):solver(s),objective(o){if(active)std::terminate();entries.clear();counters.clear();depth=0;active=true;}
 ~Root(){
  if(depth)std::terminate();
  active=false;
  std::cerr<<std::setprecision(17)<<"MATHPROFILE {\"solver\":\""<<solver<<"\",\"objective\":"<<objective<<",\"stages\":{";
  bool first=true;
  for(const auto &[n,e]:entries){std::cerr<<(first?"":",")<<"\""<<n<<"\":{\"inclusive_s\":"<<e.inclusive*1e-9<<",\"exclusive_s\":"<<e.exclusive*1e-9<<",\"calls\":"<<e.calls<<"}";first=false;}
  std::cerr<<"},\"counters\":{";first=true;
  for(const auto &[n,c]:counters){std::cerr<<(first?"":",")<<"\""<<n<<"\":"<<c;first=false;}
  std::cerr<<"}}\n";
 }
};
#endif
}
