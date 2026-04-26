#ifndef EARTHMESH_RADIAL_MESH_H
#define EARTHMESH_RADIAL_MESH_H

#include <vector>

#include <DSpecM1D/ModelConcepts>
#include <GaussQuad/All>

namespace EarthMesh {
template <typename FLOAT>
FLOAT
StandardIntervalMap(const FLOAT &x, const FLOAT &x1, const FLOAT &x2) {
   return ((x2 - x1) * x + (x1 + x2)) * 0.5;
};

// lambda to generate nodes between two radii
// template <class GQUAD>
auto
nodegenerateF(const std::vector<double> &x, double lr, double ur) {
   // vector to contain all nodes in element
   auto tmp_radius = std::vector<double>(x.size(), 0.0);
   tmp_radius[0] = lr;
   for (int idx2 = 1; idx2 < x.size() - 1; ++idx2) {
      tmp_radius[idx2] = StandardIntervalMap(x[idx2], lr, ur);
   }

   tmp_radius[x.size() - 1] = ur;
   return tmp_radius;
};

}   // namespace EarthMesh

#endif
