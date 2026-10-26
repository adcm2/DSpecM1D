#ifndef EARTHMESH_RADIALMESHDEFINITION_H
#define EARTHMESH_RADIALMESHDEFINITION_H

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace EarthMesh {

RadialMesh::RadialMesh(const double planetradius, const double ballrad,
                       const double rel_max_el, int NQ)
    : planet_radius{planetradius},
      _q{GaussQuad::GaussLobattoLegendreQuadrature1D<double>(NQ)} {
   double max_element_size = rel_max_el * planetradius;
   // lambda to generate nodes between two radii
   // auto nodegenerate = [&q = _q, NQ](double lr, double ur) {
   //    // vector to contain all nodes in element
   //    auto tmp_radius = std::vector<double>(NQ, 0.0);
   //    tmp_radius[0] = lr;
   //    for (int idx2 = 1; idx2 < NQ - 1; ++idx2) {
   //       tmp_radius[idx2] = StandardIntervalMap(q.X(idx2), lr, ur);
   //    }

   //    tmp_radius[NQ - 1] = ur;
   //    return tmp_radius;
   // };

   // vectors to store the width of each element and the number of elements per
   // layer
   std::vector<double> e_w;
   std::vector<int> N_EL;

   // internal
   N_EL.push_back(std::ceil(planet_radius / max_element_size));
   e_w.push_back(planet_radius / static_cast<double>(N_EL[0]));

   _planetary_element = std::reduce(N_EL.begin(), N_EL.end()) - 1;
   _vec_issolid.resize(_planetary_element + 1, true);

   // layer between planet and ball
   N_EL.push_back(std::ceil((ballrad - planet_radius) / max_element_size));
   e_w.push_back((ballrad - planet_radius) / static_cast<double>(N_EL[1]));

   // number of elements in total
   auto totlength = std::reduce(N_EL.begin(), N_EL.end());
   _vec_nodes.resize(totlength);
   _vec_layers.resize(totlength);

   //////////////////////////////////////////////////////////////////////
   // fill out _vec_nodes
   int idx_global = 0;
   auto x_vals = _q.Points();

   // filling planet
   for (int idxelem = 0; idxelem < N_EL[0]; ++idxelem) {
      auto lr = idxelem * e_w[0];
      auto ur = lr + e_w[0];
      if (idxelem == (N_EL[0] - 1)) {
         ur = planetradius;
      }
      _vec_nodes[idx_global] = nodegenerateF(x_vals, lr, ur);
      _vec_layers[idx_global] = 0;
      ++idx_global;
   }

   // filling layer from outer edge to ball
   for (int idxelem = 0; idxelem < N_EL[1]; ++idxelem) {
      auto lr = planetradius + idxelem * e_w[1];
      auto ur = lr + e_w[1];
      if (idxelem == (N_EL[1] - 1)) {
         ur = ballrad;
      }
      _vec_nodes[idx_global] = nodegenerateF(x_vals, lr, ur);
      _vec_layers[idx_global] = 1;
      ++idx_global;
   }
};

// setting
template <class sphericalmodel>
   requires ModelInfoConcepts::BasicSphericalDensityModel<sphericalmodel, int,
                                                          double>
RadialMesh::RadialMesh(const sphericalmodel &mymodel, int NQ,
                       const double ballrad, const double rel_max_el,
                       bool incball, double source_radius)
    : planet_radius(mymodel.OuterRadius()),
      _q{GaussQuad::GaussLobattoLegendreQuadrature1D<double>(NQ)} {
   double max_element_size = rel_max_el * mymodel.OuterRadius();
   // lambda to generate nodes between two radii
   // auto nodegenerate = [&q = _q](double lr, double ur) {
   //    // vector to contain all nodes in element
   //    auto tmp_radius = std::vector<double>(q.N(), 0.0);
   //    tmp_radius[0] = lr;
   //    for (int idx2 = 1; idx2 < q.N() - 1; ++idx2) {
   //       tmp_radius[idx2] = StandardIntervalMap(q.X(idx2), lr, ur);
   //    }
   //    tmp_radius[q.N() - 1] = ur;
   //    return tmp_radius;
   // };
   auto x_vals = _q.Points();
   // total layers and vectors to store data on elements
   int totlayers = mymodel.NumberOfLayers();
   int numlayers;
   if (incball) {
      numlayers = totlayers + 1;
   } else {
      numlayers = totlayers;
   }
   std::vector<double> e_w(numlayers);
   std::vector<int> N_EL(numlayers);

   // finding width of elements and number of elements within each layer
   for (int idx = 0; idx < totlayers; ++idx) {
      double laydepth = mymodel.UpperRadius(idx) - mymodel.LowerRadius(idx);
      N_EL[idx] = std::ceil(laydepth / max_element_size);
      e_w[idx] = laydepth / static_cast<double>(N_EL[idx]);
   };

   _planetary_element = std::reduce(N_EL.begin(), N_EL.end()) - 1;
   _vec_issolid.resize(_planetary_element + 1, true);
   // std::cout << "\nPlanetary elements: " << _planetary_element << "\n\n";
   if (incball) {
      double laydepth = ballrad - mymodel.OuterRadius();
      N_EL[totlayers] = std::ceil(laydepth / max_element_size);
      e_w[totlayers] = laydepth / static_cast<double>(N_EL[totlayers]);
   }

   // number of elements in total
   auto totlength = std::reduce(N_EL.begin(), N_EL.end());
   // std::cout << "Total elements: " << totlength << "\n\n";
   _vec_nodes.resize(totlength);
   _vec_layers.resize(totlength);

   //////////////////////////////////////////////////////////////////////
   // fill out _vec_nodes
   int idx_global = 0;
   bool currently_solid = mymodel.IsSolid(0);

   // filling up to outer edge of planet
   for (int idx = 0; idx < totlayers; ++idx) {
      // up to final element
      for (int idxelem = 0; idxelem < N_EL[idx]; ++idxelem) {
         auto lr = mymodel.LowerRadius(idx) + idxelem * e_w[idx];
         auto ur = lr + e_w[idx];
         if (idxelem == (N_EL[idx] - 1)) {
            ur = mymodel.UpperRadius(idx);
         }
         _vec_nodes[idx_global] = nodegenerateF(x_vals, lr, ur);
         _vec_layers[idx_global] = idx;
         _vec_issolid[idx_global] = mymodel.IsSolid(idx);
         if (currently_solid != mymodel.IsSolid(idx)) {
            currently_solid = mymodel.IsSolid(idx);
            _vec_fs_boundary.push_back(idx_global - 1);
         }
         if (mymodel.IsFluid(idx)) {
            _has_fluid = true;
         }

         ++idx_global;
      }
   };

   // filling layer from outer edge to ball
   if (incball) {
      // first part
      for (int idxelem = 0; idxelem < N_EL[totlayers]; ++idxelem) {
         auto lr = mymodel.OuterRadius() + idxelem * e_w[totlayers];
         auto ur = lr + e_w[totlayers];
         if (idxelem == (N_EL[totlayers] - 1)) {
            ur = ballrad;
         }
         _vec_nodes[idx_global] = nodegenerateF(x_vals, lr, ur);
         _vec_layers[idx_global] = totlayers;
         ++idx_global;
      }
   }
   if (std::isfinite(source_radius))
      splitAtRadius(source_radius);
};

inline void RadialMesh::splitAtRadius(double radius) {
   for (int e = 0; e < static_cast<int>(_vec_nodes.size()); ++e) {
      const double lo = ELR(e), hi = EUR(e);
      const double tolerance = 8.0 * std::numeric_limits<double>::epsilon() *
         std::max({1.0, std::abs(radius), std::abs(lo), std::abs(hi)});
      if (std::abs(radius - lo) <= tolerance) {
         _source_cut_radius = lo;
         return;
      }
      if (std::abs(radius - hi) <= tolerance) {
         _source_cut_radius = hi;
         return;
      }
      if (radius > lo && radius < hi) {
         auto lower = nodegenerateF(_q.Points(), lo, radius);
         auto upper = nodegenerateF(_q.Points(), radius, hi);
         const int layer = _vec_layers[e];
         const bool solid = _vec_issolid[e];
         _vec_nodes[e] = std::move(lower);
         _vec_nodes.insert(_vec_nodes.begin() + e + 1, std::move(upper));
         _vec_layers.insert(_vec_layers.begin() + e + 1, layer);
         _vec_issolid.insert(_vec_issolid.begin() + e + 1, solid);
         for (auto &boundary : _vec_fs_boundary)
            if (boundary >= e) ++boundary;
         if (e <= _planetary_element) ++_planetary_element;
         _source_cut_radius = radius;
         return;
      }
   }
   throw std::out_of_range("Source radius is outside the radial mesh.");
}
};   // namespace EarthMesh
#endif   // EARTHMESH_RADIALMESHDEFINITION_H
