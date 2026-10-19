include(FetchContent)

# Keep this dependency closure identical for in-tree and installed consumers.
set(GSHTRANS_NUMERICCONCEPTS_TAG e03a2889d43a03660d0879f46dda1406d90c9517
    CACHE STRING "Pinned NumericConcepts revision for DSpecM1D" FORCE)
set(GSHTRANS_GAUSSQUAD_TAG 7d4786d192f68077cd99a34dcdf1379ed1eeac9a
    CACHE STRING "Pinned GaussQuad revision for DSpecM1D" FORCE)
set(GSHTRANS_FFTWPP_TAG cb152c0e85f3822493c2d9e85a56307ac7b4f1c0
    CACHE STRING "Pinned FFTWpp revision for DSpecM1D" FORCE)
set(GSHTRANS_INTERPOLATION_TAG d5a1711ca6a89fc8ad8148fe13212398fe569bc0
    CACHE STRING "Pinned Interpolation revision for DSpecM1D" FORCE)

FetchContent_Declare(
    Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
    FIND_PACKAGE_ARGS 3.4.0 EXACT NAMES Eigen3)
FetchContent_Declare(
    NumericConcepts
    GIT_REPOSITORY https://github.com/da380/NumericConcepts.git
    GIT_TAG e03a2889d43a03660d0879f46dda1406d90c9517)
FetchContent_Declare(
    FFTWpp
    GIT_REPOSITORY https://github.com/da380/FFTWpp.git
    GIT_TAG cb152c0e85f3822493c2d9e85a56307ac7b4f1c0)
FetchContent_Declare(
    GaussQuad
    GIT_REPOSITORY https://github.com/da380/GaussQuad.git
    GIT_TAG 7d4786d192f68077cd99a34dcdf1379ed1eeac9a)
FetchContent_Declare(
    GSHTRANS
    GIT_REPOSITORY https://github.com/da380/GSHTrans.git
    GIT_TAG 8dd116df8df4222fe9b10a40c7213d34f2a274fb)
FetchContent_Declare(
    Interpolation
    GIT_REPOSITORY https://github.com/da380/Interpolation.git
    GIT_TAG d5a1711ca6a89fc8ad8148fe13212398fe569bc0)

function(dspecm1d_fetch_dependencies)
  set(_dspecm1d_saved_skip_install_rules ${CMAKE_SKIP_INSTALL_RULES})
  set(CMAKE_SKIP_INSTALL_RULES ON)
  FetchContent_MakeAvailable(
      Eigen3 NumericConcepts FFTWpp GaussQuad GSHTRANS Interpolation)
  set(CMAKE_SKIP_INSTALL_RULES ${_dspecm1d_saved_skip_install_rules} PARENT_SCOPE)
endfunction()
