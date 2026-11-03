#pragma once
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
namespace mathprofile {
using Clock=std::chrono::steady_clock;
struct Entry {std::uint64_t inclusive=0,exclusive=0,calls=0;};
struct Frame {const char*name=nullptr;Clock::time_point start;std::uint64_t children=0;};
inline thread_local bool active=false;
inline thread_local std::array<Frame,128> stack;
inline thread_local std::size_t depth=0;
inline thread_local std::map<std::string,Entry> entries;
inline thread_local std::map<std::string,std::uint64_t> counters;
#ifdef MATHPROFILE_DISABLED
struct Scope { explicit Scope(const char*) noexcept {} };
inline void count(const char*,std::uint64_t=1) noexcept {}
inline std::uint64_t stamp() noexcept {return 0;}
inline std::uint64_t elapsed(std::uint64_t) noexcept {return 0;}
inline void addChild(const char*,std::uint64_t,std::uint64_t=1) noexcept {}
struct Root { explicit Root(const char*,int=0) noexcept {} };
#else
struct Scope {
 bool enabled;std::size_t index=0;
 explicit Scope(const char*n):enabled(active){if(enabled){index=depth++;if(depth>stack.size())std::terminate();stack[index]={n,Clock::now(),0};}}
 ~Scope(){if(!enabled)return;const auto end=Clock::now();auto &f=stack[index];const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(end-f.start).count();--depth;auto &e=entries[f.name];e.inclusive+=elapsed;e.exclusive+=elapsed-f.children;++e.calls;if(depth)stack[depth-1].children+=elapsed;}
 Scope(const Scope&)=delete;
};
inline void count(const char*n,std::uint64_t v=1){if(active)counters[n]+=v;}
// Fine fused loops accumulate disjoint raw clock intervals locally, then
// commit once per block to avoid per-point maps, frames or new field buffers.
inline std::uint64_t stamp() noexcept {return 0;}
inline std::uint64_t elapsed(std::uint64_t) noexcept {return 0;}
inline void addChild(const char*,std::uint64_t,std::uint64_t=1) noexcept {}
struct Root {
 const char*solver;int objective;
 explicit Root(const char*s,int o=0):solver(s),objective(o){if(active)std::terminate();entries.clear();counters.clear();depth=0;active=true;}
 ~Root(){if(depth)std::terminate();active=false;std::cerr<<std::setprecision(17)<<"MATHPROFILE {\"solver\":\""<<solver<<"\",\"objective\":"<<objective<<",\"stages\":{";bool first=true;for(const auto &[n,e]:entries){std::cerr<<(first?"":",")<<"\""<<n<<"\":{\"inclusive_s\":"<<e.inclusive*1e-9<<",\"exclusive_s\":"<<e.exclusive*1e-9<<",\"calls\":"<<e.calls<<"}";first=false;}std::cerr<<"},\"counters\":{";first=true;for(const auto &[n,c]:counters){std::cerr<<(first?"":",")<<"\""<<n<<"\":"<<c;first=false;}std::cerr<<"}}\n";}
};
#endif
}
