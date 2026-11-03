#pragma once
#include <chrono>
#include <array>
#include <iostream>
namespace skprofile {
struct Entry { double seconds=0; unsigned long long calls=0; };
inline std::array<Entry,20> entries{};
inline const char* names[]={"api", "solve_modes", "prepare_workspace", "sample_fields", "pole_synthesis", "geometry_grouping", "angular_tables", "interpolation", "point_accumulation", "system_matrix", "factorization", "forward_solve", "adjoint_solve", "point_contraction", "measurement_map"};
struct Scope { int id; std::chrono::steady_clock::time_point start; Scope(int i):id(i),start(std::chrono::steady_clock::now()){} ~Scope(){entries[id].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();++entries[id].calls;} };
struct Dump { ~Dump(){std::cerr<<"SKPROFILE {"; for(int i=0;i<15;++i)std::cerr<<(i?",":"")<<"\""<<names[i]<<"\":{\"seconds\":"<<entries[i].seconds<<",\"calls\":"<<entries[i].calls<<"}"; std::cerr<<"}\n";} };
inline Dump dump;
}
