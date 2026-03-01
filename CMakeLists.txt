cmake_minimum_required (VERSION 3.14)

# set(CMAKE_C_COMPILER "gcc")
# set(CMAKE_CXX_COMPILER "g++")


project(new_coupling)
# get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_FLAGS  "-ftime-report")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2b -fconcepts-diagnostics-depth=2")
# set(CMAKE_CXX_FLAGS "-Ofast")


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/FindLAPACKE.cmake")

# Locate  FFTW3


# add in Eigen3
# add in Eigen3
# include(FetchContent)


# Helper function to fetch a dependency while disabling its install() rules.
function(FetchContent_DisableInstall TARGET_NAME)
    # In this scope, redefine the install() command to do absolutely nothing.
    macro(install)
    endmacro()

    FetchContent_MakeAvailable(${TARGET_NAME})

    # The macro override disappears automatically when this function ends.
endfunction()




include(FetchContent)
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4
)
FetchContent_MakeAvailable(Eigen3)


# find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# set(CMAKE_MODULE_PATH "/usr/local/include/eigen/cmake/;${CMAKE_MODULE_PATH}")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# Tell the Find module where to search for the SuiteSparse libraries and headers
# set(CMAKE_PREFIX_PATH "$ENV{HOME}/space/SuiteSparse" ${CMAKE_PREFIX_PATH})
# find_package(BLAS)
find_package(FFTW)
# find_package(SuiteSparse REQUIRED COMPONENTS UMFPACK)


# --- Spectra ---
FetchContent_Declare(
    Spectra
    GIT_REPOSITORY https://github.com/yixuan/spectra.git
    GIT_TAG master
)
# Call Spectra using our special function to suppress its install() command.
FetchContent_DisableInstall(Spectra)

# find_package (FFTWpp)
include(FetchContent)
FetchContent_Declare(
  FFTWpp
  GIT_REPOSITORY git@github.com:da380/FFTWpp
  GIT_TAG main
)
FetchContent_MakeAvailable(FFTWpp)

include(FetchContent)
FetchContent_Declare(
  GaussQuad
  GIT_REPOSITORY git@github.com:da380/GaussQuad
  GIT_TAG main
)
FetchContent_MakeAvailable(GaussQuad)

include(FetchContent)
FetchContent_Declare(
  PlanetaryModel
  GIT_REPOSITORY git@github.com:da380/PlanetaryModel
  GIT_TAG newdev
)
FetchContent_MakeAvailable(PlanetaryModel)

include(FetchContent)
FetchContent_Declare(
  GSHTRANS
  GIT_REPOSITORY git@github.com:da380/GSHTrans
  GIT_TAG main
)
FetchContent_MakeAvailable(GSHTRANS)


include(FetchContent)
FetchContent_Declare(
  TomographyModels
  GIT_REPOSITORY git@github.com:adcm2/TomographyModels
  GIT_TAG main
)
FetchContent_MakeAvailable(TomographyModels)

include(FetchContent)
FetchContent_Declare(
  Interpolation
  GIT_REPOSITORY git@github.com:da380/Interpolation.git
  GIT_TAG develop
)
FetchContent_MakeAvailable(Interpolation)

include(FetchContent)
FetchContent_Declare(
  gplspec
  GIT_REPOSITORY git@github.com:adcm2/gplspec.git
  GIT_TAG main
)
FetchContent_MakeAvailable(gplspec)

include(FetchContent)
FetchContent_Declare(
  EarthMesh
  GIT_REPOSITORY git@github.com:adcm2/EarthMesh.git
  GIT_TAG main
)
FetchContent_MakeAvailable(EarthMesh)

# include(FetchContent)
# FetchContent_Declare(
#   SpectraSolver
#   GIT_REPOSITORY git@github.com:adcm2/SpectraSolver.git
#   GIT_TAG develop
# )
# FetchContent_MakeAvailable(SpectraSolver)

# Set up the library
find_package(OpenMP REQUIRED)
add_library (new_coupling INTERFACE)
# target_link_libraries(new_coupling INTERFACE SuiteSparse::UMFPACK)
target_link_libraries(new_coupling INTERFACE Eigen3::Eigen)
target_link_libraries(new_coupling INTERFACE Spectra)
target_link_libraries(new_coupling INTERFACE FFTWpp)
target_link_libraries(new_coupling INTERFACE Interpolation)
target_link_libraries(new_coupling INTERFACE GaussQuad)
target_link_libraries(new_coupling INTERFACE PlanetaryModel)
target_link_libraries(new_coupling INTERFACE TomographyModels)
target_link_libraries(new_coupling INTERFACE gplspec)
target_link_libraries(new_coupling INTERFACE EarthMesh)
target_link_libraries(new_coupling INTERFACE OpenMP::OpenMP_CXX)
# target_link_libraries(new_coupling INTERFACE SpectraSolver)
# target_link_libraries(new_coupling INTERFACE GSHTRANS)
target_include_directories (new_coupling INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
)
# set_target_properties(new_coupling PROPERTIES UNITY_BUILD ON)
# target_precompile_headers(new_coupling INTERFACE <vector> <string> <iostream> "pch.h")

# blas and lapack
# find_package(OpenMP REQUIRED)





# optionally add in the examples and tests
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    # We're in the root, define additional targets for developers.
    option(MY_PROJECT_BUILD_EXAMPLES   "whether or not examples should be built" ON)
    # option(MY_PROJECT_BUILD_TESTS      "whether or not tests should be built" ON)


    
    # include(FetchContent)
    # FetchContent_Declare(
    #   Interpolation
    #   GIT_REPOSITORY https://github.com/da380/Interpolation.git
    #   GIT_TAG main
    #   )
    # FetchContent_MakeAvailable(Interpolation)

    if(MY_PROJECT_BUILD_EXAMPLES)
        # add_subdirectory(examples)
        # add_subdirectory(spheroidal_examples)
        add_subdirectory(benchmarks)
    endif()
    # if(MY_PROJECT_BUILD_TESTS)
    #     enable_testing()
    #     add_subdirectory(tests)
    # endif()
endif()

# add_subdirectory(Gravitational_Field)
# add_subdirectory(examples)
