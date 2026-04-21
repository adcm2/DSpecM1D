cmake_minimum_required(VERSION 3.24)
project(DSpecM1D_Draft VERSION 0.1.0 LANGUAGES CXX)

# --- Core Settings ---
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new build directory.")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# --- System Dependencies ---
find_package(FFTW)
find_package(OpenMP REQUIRED)

# --- Parallel Dependency Fetching (CMake 3.24+) ---
# Store downloads in a hidden .deps folder so 'rm -rf build' doesn't delete them!
set(FETCHCONTENT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.deps)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL "Disable FetchContent updates")
include(FetchContent)

# 1. Declare everything upfront.
# The FIND_PACKAGE_ARGS tells CMake to search your system for Eigen3 first.
# If it finds it locally, it SKIPS the git download entirely!
FetchContent_Declare(
    Eigen3           
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git         
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
    FIND_PACKAGE_ARGS NAMES Eigen3
)
FetchContent_Declare(FFTWpp           GIT_REPOSITORY https://github.com/da380/FFTWpp.git           GIT_TAG main   GIT_SHALLOW TRUE)
FetchContent_Declare(GaussQuad        GIT_REPOSITORY https://github.com/da380/GaussQuad.git        GIT_TAG main   GIT_SHALLOW TRUE)
FetchContent_Declare(PlanetaryModel   GIT_REPOSITORY https://github.com/da380/PlanetaryModel.git   GIT_TAG newdev GIT_SHALLOW TRUE)
FetchContent_Declare(GSHTRANS         GIT_REPOSITORY https://github.com/da380/GSHTrans.git         GIT_TAG main   GIT_SHALLOW TRUE)
FetchContent_Declare(TomographyModels GIT_REPOSITORY https://github.com/adcm2/TomographyModels.git GIT_TAG main   GIT_SHALLOW TRUE)
FetchContent_Declare(Interpolation    GIT_REPOSITORY https://github.com/da380/Interpolation.git    GIT_TAG develop GIT_SHALLOW TRUE)
FetchContent_Declare(EarthMesh        GIT_REPOSITORY https://github.com/adcm2/EarthMesh.git        GIT_TAG main   GIT_SHALLOW TRUE)
FetchContent_Declare(SpectraSolver    GIT_REPOSITORY https://github.com/adcm2/SpectraSolver.git    GIT_TAG develop GIT_SHALLOW TRUE)

# 2. Modify our helper function to accept a LIST of targets
function(FetchContent_MakeAvailable_NoInstall)
    macro(install)
    endmacro()
    FetchContent_MakeAvailable(${ARGV})
endfunction()

# 3. Fire them all off at once! 
# CMake will handle the "local vs git" logic automatically.
FetchContent_MakeAvailable_NoInstall(
    Eigen3 FFTWpp GaussQuad PlanetaryModel GSHTRANS 
    TomographyModels Interpolation EarthMesh SpectraSolver
)

# --- Target Definition ---
add_library(DSpecM1D INTERFACE)

target_link_libraries(DSpecM1D INTERFACE 
    Eigen3::Eigen
    FFTWpp
    Interpolation
    GaussQuad
    PlanetaryModel
    TomographyModels
    EarthMesh
    SpectraSolver
    OpenMP::OpenMP_CXX
)

target_include_directories(DSpecM1D INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# --- Benchmarks ---
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    enable_testing()
    option(DSPECM1D_BUILD_BENCHMARKS "Whether or not benchmarks should be built" ON)
    option(DSPECM1D_BUILD_TUTORIALS "Whether or not tutorials should be built" ON)
    option(DSPECM1D_BUILD_TESTS "Build DSPECM1D Tests" OFF)
    option(DSPECM1D_ENABLE_SMOKE_TESTS "Register lightweight executable smoke tests" OFF)

    if(DSPECM1D_BUILD_BENCHMARKS OR DSPECM1D_BUILD_TUTORIALS)
        add_custom_target(CopyData
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/data
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )
    endif()

    if(DSPECM1D_BUILD_BENCHMARKS)
        add_subdirectory(benchmarks)
    endif()

    if(DSPECM1D_BUILD_TUTORIALS)
        add_subdirectory(tutorials)
    endif()

    if(DSPECM1D_BUILD_TESTS)
        add_subdirectory(tests)
    endif()
endif()