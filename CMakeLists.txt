cmake_minimum_required(VERSION 3.24)
project(DSpecM1D_Draft VERSION 0.1.0 LANGUAGES CXX)

# --- Core Settings ---
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new build directory.")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# --- System Dependencies ---
find_package(OpenMP REQUIRED)

# --- Parallel Dependency Fetching (CMake 3.24+) ---
# Keep FetchContent build/subbuild state local to each build tree. Sharing this
# cache across build directories can mix generators and break reconfigure.
set(FETCHCONTENT_BASE_DIR "${CMAKE_BINARY_DIR}/_deps" CACHE PATH
    "FetchContent dependency cache for this build tree" FORCE)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL "Disable FetchContent updates")
include(FetchContent)
include(cmake/DSpecM1DDependencies.cmake)
dspecm1d_fetch_dependencies()

# --- Target Definition ---
add_library(DSpecM1D INTERFACE)

target_link_libraries(DSpecM1D INTERFACE 
    Eigen3::Eigen
    NumericConcepts
    FFTWpp
    GaussQuad
    Interpolation
    GSHTrans
    OpenMP::OpenMP_CXX
)

target_compile_features(DSpecM1D INTERFACE cxx_std_23)

target_include_directories(DSpecM1D INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/DSpecM1D
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DSpecM1DDependencies.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DSpecM1D
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/DSpecM1DConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DSpecM1DConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DSpecM1DConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DSpecM1D
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/DSpecM1DConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/DSpecM1DConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DSpecM1D
)

# --- Top-Level Build Options ---
option(DSPECM1D_BUILD_BENCHMARKS "Whether or not benchmarks should be built" ON)
option(DSPECM1D_BUILD_TUTORIALS "Whether or not tutorials should be built" ON)
option(DSPECM1D_BUILD_TOOLS "Build DSpecM1D command-line tools" ON)
option(DSPECM1D_BUILD_TESTS "Build DSpecM1D tests" OFF)
option(DSPECM1D_BUILD_LOVE_NUMBERS "Build the static Love-number module" OFF)
option(DSPECM1D_ENABLE_SMOKE_TESTS "Register lightweight executable smoke tests" OFF)
option(DSPECM1D_ENABLE_MIGRATION_TESTS "Register temporary migration validation checks" OFF)
option(DSPECM1D_ENABLE_PAPER_VALIDATION "Register optional paper reproduction validation checks" OFF)
option(DSPECM1D_BUILD_DOCS "Build DSpecM1D documentation targets" OFF)

# --- Developer Targets ---
if(PROJECT_IS_TOP_LEVEL)
    enable_testing()

    if(DSPECM1D_BUILD_BENCHMARKS OR DSPECM1D_BUILD_TUTORIALS OR DSPECM1D_BUILD_TOOLS)
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

    if(DSPECM1D_BUILD_TOOLS)
        add_subdirectory(tools)
    endif()

    if(DSPECM1D_BUILD_TESTS)
        add_subdirectory(tests)
    endif()

    if(DSPECM1D_BUILD_DOCS)
        add_subdirectory(docs)
    endif()
endif()

if(DSPECM1D_BUILD_LOVE_NUMBERS)
    add_subdirectory(love_numbers)
endif()
