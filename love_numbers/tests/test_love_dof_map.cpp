#include <gtest/gtest.h>

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/SEM/SEM.h>

#include "LoveDofMap.h"

namespace {

using DSpecM1D::LoveNumbers::detail::Field;
using DSpecM1D::LoveNumbers::detail::LoveDofMap;

EarthMesh::RadialMesh makeAllSolidMesh() {
  return EarthMesh::RadialMesh(1.0, 1.0, 0.25, 4);
}

EarthMesh::RadialMesh makeInternalFluidMesh() {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL};
  return EarthMesh::RadialMesh(model, 4, 1.0, 0.2, false);
}

int uniqueRadialNodes(const EarthMesh::RadialMesh &mesh) {
  return mesh.NE() * (mesh.NN() - 1) + 1;
}

int uniqueSolidNodes(const EarthMesh::RadialMesh &mesh) {
  int nodes = 0;
  for (int element = 0; element < mesh.NE(); ++element) {
    if (!mesh.IsSolid(element)) {
      continue;
    }
    nodes += mesh.NN() - 1;
    if (element == 0 || !mesh.IsSolid(element - 1)) {
      ++nodes;
    }
  }
  return nodes;
}

TEST(LoveDofMapTests, DegreeZeroInterleavesEveryNodeAndMatchesRadialSemMap) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL};
  const Full1D::SEM sem(model, 0.2, 4, 2);
  const auto &mesh = sem.mesh();
  const LoveDofMap map(mesh, 0);

  EXPECT_EQ(map.size(), 2 * uniqueRadialNodes(mesh));
  for (int element = 0; element < mesh.NE(); ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      const int radial_node = element * (mesh.NN() - 1) + node;
      ASSERT_TRUE(map.index(Field::U, element, node));
      ASSERT_TRUE(map.index(Field::P, element, node));
      EXPECT_FALSE(map.index(Field::V, element, node));
      EXPECT_EQ(*map.index(Field::U, element, node), 2 * radial_node);
      EXPECT_EQ(*map.index(Field::P, element, node), 2 * radial_node + 1);
      EXPECT_EQ(*map.index(Field::U, element, node),
                sem.ltgR(0, element, node));
      EXPECT_EQ(*map.index(Field::P, element, node),
                sem.ltgR(1, element, node));
    }
  }
}

TEST(LoveDofMapTests, DegreeZeroRetainsCentreFields) {
  const auto mesh = makeInternalFluidMesh();
  const LoveDofMap map(mesh, 0);

  ASSERT_TRUE(map.index(Field::U, 0, 0));
  ASSERT_TRUE(map.index(Field::P, 0, 0));
  EXPECT_EQ(*map.index(Field::U, 0, 0), 0);
  EXPECT_EQ(*map.index(Field::P, 0, 0), 1);
}

TEST(LoveDofMapTests, AllSolidDegreeTwoInterleavesEveryField) {
  const auto mesh = makeAllSolidMesh();
  const LoveDofMap map(mesh, 2);

  EXPECT_EQ(map.size(), 3 * uniqueRadialNodes(mesh));
  for (int element = 0; element < mesh.NE(); ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      const int radial_node = element * (mesh.NN() - 1) + node;
      EXPECT_EQ(*map.index(Field::U, element, node), 3 * radial_node);
      EXPECT_EQ(*map.index(Field::V, element, node), 3 * radial_node + 1);
      EXPECT_EQ(*map.index(Field::P, element, node), 3 * radial_node + 2);
    }
  }
}

TEST(LoveDofMapTests, SolidCentreRetainsThreeIndependentFields) {
  const auto mesh = makeAllSolidMesh();
  const LoveDofMap map(mesh, 2);

  const auto u = map.index(Field::U, 0, 0);
  const auto v = map.index(Field::V, 0, 0);
  const auto p = map.index(Field::P, 0, 0);
  ASSERT_TRUE(u);
  ASSERT_TRUE(v);
  ASSERT_TRUE(p);
  EXPECT_NE(*u, *v);
  EXPECT_NE(*u, *p);
  EXPECT_NE(*v, *p);
}

TEST(LoveDofMapTests, PotentialIsSharedAcrossEveryElementBoundary) {
  const auto mesh = makeInternalFluidMesh();
  const LoveDofMap map(mesh, 2);

  for (int element = 0; element + 1 < mesh.NE(); ++element) {
    EXPECT_EQ(map.index(Field::P, element, mesh.NN() - 1),
              map.index(Field::P, element + 1, 0));
  }
}

TEST(LoveDofMapTests, DisplacementsAreSharedAcrossSolidBoundaries) {
  const auto mesh = makeInternalFluidMesh();
  const LoveDofMap map(mesh, 2);
  bool checked_boundary = false;

  for (int element = 0; element + 1 < mesh.NE(); ++element) {
    if (!mesh.IsSolid(element) || !mesh.IsSolid(element + 1)) {
      continue;
    }
    EXPECT_EQ(map.index(Field::U, element, mesh.NN() - 1),
              map.index(Field::U, element + 1, 0));
    EXPECT_EQ(map.index(Field::V, element, mesh.NN() - 1),
              map.index(Field::V, element + 1, 0));
    checked_boundary = true;
  }

  EXPECT_TRUE(checked_boundary);
}

TEST(LoveDofMapTests, DisplacementsAreAbsentFromFluidElements) {
  const auto mesh = makeInternalFluidMesh();
  const LoveDofMap map(mesh, 2);
  bool checked_fluid = false;

  for (int element = 0; element < mesh.NE(); ++element) {
    if (!mesh.IsFluid(element)) {
      continue;
    }
    for (int node = 0; node < mesh.NN(); ++node) {
      EXPECT_FALSE(map.index(Field::U, element, node));
      EXPECT_FALSE(map.index(Field::V, element, node));
      EXPECT_TRUE(map.index(Field::P, element, node));
    }
    checked_fluid = true;
  }

  EXPECT_TRUE(checked_fluid);
}

TEST(LoveDofMapTests, InterfacesKeepDisplacementsOnlyOnSolidSide) {
  const auto mesh = makeInternalFluidMesh();
  const LoveDofMap map(mesh, 2);
  bool checked_solid_to_fluid = false;
  bool checked_fluid_to_solid = false;

  for (int lower = 0; lower + 1 < mesh.NE(); ++lower) {
    const int upper = lower + 1;
    if (mesh.IsSolid(lower) && mesh.IsFluid(upper)) {
      EXPECT_TRUE(map.index(Field::U, lower, mesh.NN() - 1));
      EXPECT_TRUE(map.index(Field::V, lower, mesh.NN() - 1));
      EXPECT_FALSE(map.index(Field::U, upper, 0));
      EXPECT_FALSE(map.index(Field::V, upper, 0));
      checked_solid_to_fluid = true;
    }
    if (mesh.IsFluid(lower) && mesh.IsSolid(upper)) {
      EXPECT_FALSE(map.index(Field::U, lower, mesh.NN() - 1));
      EXPECT_FALSE(map.index(Field::V, lower, mesh.NN() - 1));
      EXPECT_TRUE(map.index(Field::U, upper, 0));
      EXPECT_TRUE(map.index(Field::V, upper, 0));
      checked_fluid_to_solid = true;
    }
  }

  EXPECT_TRUE(checked_solid_to_fluid);
  EXPECT_TRUE(checked_fluid_to_solid);
}

TEST(LoveDofMapTests, InternalFluidDimensionsMatchTopologyFormulas) {
  const auto mesh = makeInternalFluidMesh();
  const int radial_nodes = uniqueRadialNodes(mesh);
  const int solid_nodes = uniqueSolidNodes(mesh);

  EXPECT_EQ(LoveDofMap(mesh, 0).size(), 2 * radial_nodes);
  EXPECT_EQ(LoveDofMap(mesh, 1).size(),
            radial_nodes + 2 * solid_nodes - 1);
  EXPECT_EQ(LoveDofMap(mesh, 2).size(),
            radial_nodes + 2 * solid_nodes);
}

TEST(LoveDofMapTests, DegreeOneOnlyOmitsSurfacePotential) {
  const auto mesh = makeInternalFluidMesh();
  const LoveDofMap degree_one(mesh, 1);
  const LoveDofMap degree_two(mesh, 2);
  const int surface_element = mesh.NE() - 1;
  const int surface_node = mesh.NN() - 1;

  EXPECT_EQ(degree_one.size(), degree_two.size() - 1);
  for (int element = 0; element < mesh.NE(); ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      EXPECT_EQ(degree_one.index(Field::U, element, node),
                degree_two.index(Field::U, element, node));
      EXPECT_EQ(degree_one.index(Field::V, element, node),
                degree_two.index(Field::V, element, node));
      if (element != surface_element || node != surface_node) {
        EXPECT_EQ(degree_one.index(Field::P, element, node),
                  degree_two.index(Field::P, element, node));
      }
    }
  }
  EXPECT_FALSE(degree_one.index(Field::P, surface_element, surface_node));
  EXPECT_TRUE(degree_two.index(Field::P, surface_element, surface_node));
}

TEST(LoveDofMapTests, DegreeOneRetainsSurfaceDisplacements) {
  const auto mesh = makeInternalFluidMesh();
  const LoveDofMap map(mesh, 1);
  const int surface_element = mesh.NE() - 1;
  const int surface_node = mesh.NN() - 1;

  EXPECT_TRUE(map.index(Field::U, surface_element, surface_node));
  EXPECT_TRUE(map.index(Field::V, surface_element, surface_node));
  EXPECT_FALSE(map.index(Field::P, surface_element, surface_node));
}

TEST(LoveDofMapTests, RejectsNegativeDegree) {
  const auto mesh = makeInternalFluidMesh();

  EXPECT_THROW(LoveDofMap(mesh, -1), std::invalid_argument);
}

}   // namespace
