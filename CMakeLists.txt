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
include(FetchContent)

# 1. Declare everything upfront. (No downloading happens here).
FetchContent_Declare(Eigen3           GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git         GIT_TAG 3.4.0  GIT_SHALLOW TRUE)
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
    # ${ARGV} expands to all arguments passed to this function
    FetchContent_MakeAvailable(${ARGV})
endfunction()

# 3. Fire them all off at once! CMake will download these in parallel.
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
    $<INSTALL_INTERFACE:include>
)

# --- Install Rules (header-only INTERFACE target) ---
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Export the interface target
install(TARGETS DSpecM1D
    EXPORT DSpecM1DTargets
)

# Install all public headers under include/DSpecM1D/ so that
# #include <DSpecM1D/All> resolves correctly after installation.
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/DSpecM1D/src/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/DSpecM1D
    FILES_MATCHING PATTERN "*.h"
)

# Install umbrella headers (extensionless)
install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/DSpecM1D/All
    ${CMAKE_CURRENT_SOURCE_DIR}/DSpecM1D/Timer
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/DSpecM1D
)

# Install implementation/public headers preserving src/ path
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/DSpecM1D/src/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/DSpecM1D/src
    FILES_MATCHING PATTERN "*.h"
)

# Write the targets file
install(EXPORT DSpecM1DTargets
    FILE      DSpecM1DTargets.cmake
    NAMESPACE DSpecM1D::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DSpecM1D
)

# Version file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/DSpecM1DConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Config file from template
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DSpecM1DConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/DSpecM1DConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DSpecM1D
    NO_SET_AND_CHECK_MACRO   # not needed for header-only; avoids unused-macro warnings
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/DSpecM1DConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/DSpecM1DConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DSpecM1D
)

# --- Benchmarks ---
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(DSPECM1D_BUILD_BENCHMARKS "Whether or not benchmarks should be built" ON)
    option(DSPECM1D_BUILD_TUTORIALS "Whether or not tutorials should be built" ON)

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
endif()

# --- Tests (must come after benchmarks so ex3 target exists) ---
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    enable_testing()
    option(DSPECM1D_BUILD_TESTS "Build CTest smoke tests" ON)
    if(DSPECM1D_BUILD_TESTS)
        add_subdirectory(tests)
    endif()
endif()