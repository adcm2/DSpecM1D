#include <gtest/gtest.h>
#include <PlanetaryModel/All>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SourceInfo.h>
#include "test_utils.h"

TEST(SEMComponentTests, LocalToGlobalMapsAreMonotonicOnMinimalModel) {
  prem_norm<double> norm;
  auto model =
      EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm, "true");
  Full1D::SEM sem(model, 0.05, 4, 4);

  EXPECT_LT(sem.ltgS(0, 0, 0), sem.ltgS(1, 0, 0));
  EXPECT_LT(sem.ltgS(1, 0, 0), sem.ltgS(2, 0, 0));
  EXPECT_LT(sem.ltgR(0, 0, 0), sem.ltgR(1, 0, 0));
  EXPECT_LE(sem.el(), sem.eu());
}
