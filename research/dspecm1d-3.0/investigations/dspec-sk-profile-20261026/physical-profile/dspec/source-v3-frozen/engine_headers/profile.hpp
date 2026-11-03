#pragma once
#include <chrono>
#include <map>
#include <string>
#include <iostream>
namespace auditprofile {
#ifdef DSFAST_DISABLE_PROFILE
struct Entry {double seconds=0; unsigned long long calls=0;};
struct NullValue {
  double seconds=0; unsigned long long calls=0;
  NullValue &operator++() noexcept { return *this; }
  template<class T> NullValue &operator+=(const T&) noexcept { return *this; }
};
struct NullMap { NullValue &operator[](const char*) { static NullValue value; return value; } };
inline NullMap data,counters;
struct Scope { explicit Scope(const char*) noexcept {} };
struct Dump { };
inline Dump dump;
#else
struct Entry {double seconds=0; unsigned long long calls=0;};
inline std::map<std::string,Entry> data;
inline std::map<std::string,unsigned long long> counters;
struct Scope {const char* name;std::chrono::steady_clock::time_point start; Scope(const char* n):name(n),start(std::chrono::steady_clock::now()){} ~Scope(){auto& e=data[name];e.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();++e.calls;}};
struct Dump {~Dump(){std::cerr<<"AUDITPROFILE {\"timers\":{";bool first=true;for(auto&[n,e]:data){std::cerr<<(first?"":",")<<"\""<<n<<"\":{\"seconds\":"<<e.seconds<<",\"calls\":"<<e.calls<<"}";first=false;}std::cerr<<"},\"counters\":{";first=true;for(auto&[n,v]:counters){std::cerr<<(first?"":",")<<"\""<<n<<"\":"<<v;first=false;}std::cerr<<"}}\n";}};inline Dump dump;
#endif
}
