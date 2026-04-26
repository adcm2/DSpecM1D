#ifndef EARTHMESH_RADIALMESHDECLARATION_H
#define EARTHMESH_RADIALMESHDECLARATION_H

namespace EarthMesh {
class RadialMesh {
 public:
   // constructors
   RadialMesh() {};   // default

   // simple constructor
   RadialMesh(const double, const double, const double, int);

   // constructor for spherical density model
   //

   template <class sphericalmodel>
      requires ModelInfoConcepts::BasicSphericalDensityModel<sphericalmodel,
                                                             int, double>
   RadialMesh(const sphericalmodel &, int, const double, const double,
              bool = true);

   // return functions of nodes
   auto NE() const { return _vec_nodes.size(); };   // #elements
   auto NN() const { return _q.N(); };              // #nodes per element
   auto NodeRadius(int idxelem, int idxnode) const {
      return _vec_nodes[idxelem][idxnode];
   };

   // information on planet
   auto LayerNumber(int idxelem) const { return _vec_layers[idxelem]; };
   auto NL() const { return _vec_layers.back() + 1; };   // #layers

   // element lower radius
   auto ELR(int idxelem) const { return _vec_nodes[idxelem][0]; };

   // element upper radius
   auto EUR(int idxelem) const { return _vec_nodes[idxelem].back(); };

   // element width
   auto EW(int idxelem) const {
      return _vec_nodes[idxelem].back() - _vec_nodes[idxelem][0];
   };

   // planet radius
   auto PR() const { return planet_radius; };

   // outer radius
   auto OR() const { return _vec_nodes.back().back(); };

   // fluid or not:
   auto IsSolid(int idxelem) const {
      if (idxelem < _vec_issolid.size()) {
         return _vec_issolid[idxelem];
      } else {
         return false;
      }
   };
   auto IsFluid(int idxelem) const {
      if (idxelem < _vec_issolid.size()) {
         return !_vec_issolid[idxelem];
      } else {
         return false;
      }
   };
   auto HasFluid() const { return _has_fluid; };

   // fluid-solid boundaries
   auto FS_Boundaries() const { return _vec_fs_boundary; };
   const std::vector<int> &FS_BoundariesP() const { return _vec_fs_boundary; };

   // outer planetary element
   auto OuterPlanetaryElement() const { return _planetary_element; };

   // mesh
   GaussQuad::Quadrature1D<double> &GLL() { return _q; };

 private:
   std::vector<int> _vec_layers;
   std::vector<std::vector<double>> _vec_nodes;
   double planet_radius;
   GaussQuad::Quadrature1D<double> _q;
   std::vector<bool> _vec_issolid;
   bool _has_fluid = false;
   int _planetary_element;
   std::vector<int> _vec_fs_boundary;
};
};   // namespace EarthMesh
#endif   // EARTHMESH_RADIALMESHDEFINITION_H
