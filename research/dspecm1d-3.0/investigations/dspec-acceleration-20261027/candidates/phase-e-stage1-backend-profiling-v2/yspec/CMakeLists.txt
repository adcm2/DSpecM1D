# 3.24 for FetchContent_Declare(... FIND_PACKAGE_ARGS), with which every
# dependency below is declared.
cmake_minimum_required(VERSION 3.24)

project(YSpec
        VERSION 0.1.0
        DESCRIPTION "Synthetic seismograms for spherically symmetric earth models by direct radial integration"
        LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not supported; use cmake -B build.")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(FetchContent)

# Dependencies.
#
# All are David Al-Attar's header-only libraries, looked for on the system
# first and fetched otherwise. They are declared here, before any of them is
# made available, so that these pins win over the ones the libraries declare
# for each other: FetchContent keeps the first declaration of a name. Each pin
# is a cache variable, so a build can move one without editing this file, and
# a local checkout can be used instead with
#
#   -DFETCHCONTENT_SOURCE_DIR_<NAME>=/path/to/checkout
#
# (see CMakePresets.json for a preset that does this for ~/dev).
set(YSPEC_NUMERICCONCEPTS_TAG v1.0.0 CACHE STRING "NumericConcepts revision")
set(YSPEC_INTERPOLATION_TAG v1.1.2 CACHE STRING "Interpolation revision")
set(YSPEC_FFTWPP_TAG v1.0.1 CACHE STRING "FFTWpp revision")
set(YSPEC_GAUSSQUAD_TAG 7d4786d192f68077cd99a34dcdf1379ed1eeac9a
    CACHE STRING "GaussQuad revision")
set(YSPEC_GSHTRANS_TAG 8dd116df8df4222fe9b10a40c7213d34f2a274fb
    CACHE STRING "GSHTrans revision")
set(YSPEC_PLANETMODEL_TAG v1.3.0 CACHE STRING "planetmodel revision")
set(YSPEC_POCKETFFT_TAG 4d865cba24359b05b37da4e75aee24498c3b78e7
    CACHE STRING "pocketfft revision")
mark_as_advanced(YSPEC_NUMERICCONCEPTS_TAG YSPEC_INTERPOLATION_TAG
                 YSPEC_FFTWPP_TAG YSPEC_GAUSSQUAD_TAG YSPEC_GSHTRANS_TAG
                 YSPEC_PLANETMODEL_TAG YSPEC_POCKETFFT_TAG)

# The one Fourier transform YSpec makes, of each spectrum to its time series:
# by FFTW, through FFTWpp, or by pocketfft (header-only, the transform inside
# numpy). FFTW is GPL and pocketfft BSD, so the Python wheels, which carry the
# library they use, take pocketfft and stay BSD. Nothing else uses FFTW:
# GSHTrans gives only its Wigner functions, which need none, and without FFTW
# its headers are used directly rather than through its own CMake, which
# would look for FFTW.
set(YSPEC_FFT FFTW CACHE STRING "The FFT library: FFTW or pocketfft")
set_property(CACHE YSPEC_FFT PROPERTY STRINGS FFTW pocketfft)
if(NOT YSPEC_FFT MATCHES "^(FFTW|pocketfft)$")
  message(FATAL_ERROR "YSPEC_FFT is FFTW or pocketfft, not ${YSPEC_FFT}")
endif()

FetchContent_Declare(NumericConcepts
  GIT_REPOSITORY https://github.com/da380/NumericConcepts.git
  GIT_TAG ${YSPEC_NUMERICCONCEPTS_TAG}
  FIND_PACKAGE_ARGS)
FetchContent_Declare(Interpolation
  GIT_REPOSITORY https://github.com/da380/Interpolation.git
  GIT_TAG ${YSPEC_INTERPOLATION_TAG}
  FIND_PACKAGE_ARGS)
FetchContent_Declare(GaussQuad
  GIT_REPOSITORY https://github.com/da380/GaussQuad.git
  GIT_TAG ${YSPEC_GAUSSQUAD_TAG}
  FIND_PACKAGE_ARGS)
if(YSPEC_FFT STREQUAL "FFTW")
  FetchContent_Declare(FFTWpp
    GIT_REPOSITORY https://github.com/da380/FFTWpp.git
    GIT_TAG ${YSPEC_FFTWPP_TAG}
    FIND_PACKAGE_ARGS)
  FetchContent_Declare(GSHTrans
    GIT_REPOSITORY https://github.com/da380/GSHTrans.git
    GIT_TAG ${YSPEC_GSHTRANS_TAG}
    FIND_PACKAGE_ARGS)
else()
  # Fetched but not added: a SOURCE_SUBDIR without a CMakeLists.txt.
  FetchContent_Declare(GSHTrans
    GIT_REPOSITORY https://github.com/da380/GSHTrans.git
    GIT_TAG ${YSPEC_GSHTRANS_TAG}
    SOURCE_SUBDIR headers-only)
  FetchContent_Declare(pocketfft
    GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
    GIT_TAG ${YSPEC_POCKETFFT_TAG}
    SOURCE_SUBDIR headers-only)
endif()
# The C++ library is a subdirectory of a repository that is mostly Python.
FetchContent_Declare(planetmodel
  GIT_REPOSITORY https://github.com/da380/planetmodel.git
  GIT_TAG ${YSPEC_PLANETMODEL_TAG}
  SOURCE_SUBDIR cpp)

# Three-dimensional kernels are written as planetmodel's model files, over
# netCDF-C; without it YSpec does everything else. GSHTrans's matrix kernel,
# the only part needing BLAS, is not used.
option(YSPEC_WITH_NETCDF
       "Write three-dimensional kernels as model files (needs netCDF-C)" ON)
set(PLANETMODEL_WITH_NETCDF ${YSPEC_WITH_NETCDF} CACHE BOOL "" FORCE)
set(GSHTRANS_WITH_BLAS OFF CACHE STRING "")

if(YSPEC_FFT STREQUAL "FFTW")
  FetchContent_MakeAvailable(NumericConcepts Interpolation FFTWpp GaussQuad
                             GSHTrans planetmodel)
else()
  FetchContent_MakeAvailable(NumericConcepts Interpolation GaussQuad GSHTrans
                             pocketfft planetmodel)
  add_library(YSpecGSHTrans INTERFACE)
  target_include_directories(YSpecGSHTrans SYSTEM INTERFACE
    ${gshtrans_SOURCE_DIR})
  target_link_libraries(YSpecGSHTrans INTERFACE
    NumericConcepts::NumericConcepts)
  add_library(GSHTrans::GSHTrans ALIAS YSpecGSHTrans)
  add_library(YSpecPocketfft INTERFACE)
  target_include_directories(YSpecPocketfft SYSTEM INTERFACE
    ${pocketfft_SOURCE_DIR})
  # YSpec transforms on its own threads.
  target_compile_definitions(YSpecPocketfft INTERFACE
    POCKETFFT_NO_MULTITHREADING YSPEC_WITH_POCKETFFT)
endif()

# Every (l, omega) solve is independent, and the loop over them is threaded.
# Off, the same code runs serially and gives the same answers.
option(YSPEC_WITH_OPENMP "Thread with OpenMP" ON)
if(YSPEC_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

option(YSPEC_WARNINGS_AS_ERRORS "Treat warnings in our own code as errors" OFF)
set(YSPEC_WARNINGS -Wall -Wextra)
if(YSPEC_WARNINGS_AS_ERRORS)
  list(APPEND YSPEC_WARNINGS -Werror)
endif()

# What the library is built from but does not show: the numerics and the
# libraries they use. The tests of the internals link this too.
add_library(YSpecInternal INTERFACE)
target_include_directories(YSpecInternal INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(YSpecInternal INTERFACE
  NumericConcepts::NumericConcepts
  Interpolation::Interpolation
  GaussQuad::GaussQuad
  GSHTrans::GSHTrans
  planetmodel::planetmodel)
if(YSPEC_FFT STREQUAL "FFTW")
  target_link_libraries(YSpecInternal INTERFACE FFTWpp::FFTWpp)
else()
  target_link_libraries(YSpecInternal INTERFACE YSpecPocketfft)
endif()

# The library. Its public headers, in YSpec/, use nothing but the standard
# library; everything else is private. Static unless BUILD_SHARED_LIBS is on,
# and position-independent either way so that it can go into a Python module.
add_library(YSpec
  src/AngularFunctions.cpp
  src/Coefficients.cpp
  src/Input.cpp
  src/InputSections.cpp
  src/Instrument.cpp
  src/Kernels.cpp
  src/Location.cpp
  src/Model.cpp
  src/PhaseEProfile.cpp
  src/RadialSolver.cpp
  src/Run.cpp
  src/Snapshots.cpp
  src/Synthetics.cpp
  src/TimeSeries.cpp
  src/SampledKernels.cpp)
add_library(YSpec::YSpec ALIAS YSpec)
set_target_properties(YSpec PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(YSpec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_compile_features(YSpec PUBLIC cxx_std_23)
target_link_libraries(YSpec PRIVATE YSpecInternal)
target_compile_options(YSpec PRIVATE ${YSPEC_WARNINGS})
# The version, which the files the library writes record.
target_compile_definitions(YSpec PRIVATE YSPEC_VERSION="${PROJECT_VERSION}")
if(YSPEC_WITH_OPENMP)
  target_link_libraries(YSpec PRIVATE OpenMP::OpenMP_CXX)
endif()

# The programs, one for each kind of output. A Python wheel has its own,
# in Python, and builds none of these.
option(YSPEC_BUILD_PROGRAMS "Build the programs" ON)
if(YSPEC_BUILD_PROGRAMS)
  set(YSPEC_PROGRAMS yspec yspec-snapshot yspec-kernels)
  foreach(program IN LISTS YSPEC_PROGRAMS)
    add_executable(${program} app/${program}.cpp)
    target_link_libraries(${program} PRIVATE YSpec::YSpec)
    target_compile_definitions(${program} PRIVATE
      YSPEC_VERSION="${PROJECT_VERSION}")
    target_compile_options(${program} PRIVATE ${YSPEC_WARNINGS})
    if(YSPEC_WITH_OPENMP)
      target_link_libraries(${program} PRIVATE OpenMP::OpenMP_CXX)
    endif()
  endforeach()

  # Only the programs are installed for now; exporting the library as a CMake
  # package has to wait until planetmodel exports one of its own.
  include(GNUInstallDirs)
  install(TARGETS ${YSPEC_PROGRAMS}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# The Python module, yspec._core, which the package in python/ wraps. The
# wheels are built through pyproject.toml, which turns this on.
option(YSPEC_PYTHON "Build the Python module" OFF)
if(YSPEC_PYTHON)
  add_subdirectory(python)
endif()

# The examples -- inputs, models and the scripts that run them -- copied into
# the build tree, so that from the build directory
# `bin/yspec examples/input.json`, or `examples/record-section/run.sh`, runs
# and writes its output there rather than in the sources. They are copied
# again when they change.
option(YSPEC_BUILD_EXAMPLES "Copy the examples into the build tree" ON)
if(PROJECT_IS_TOP_LEVEL AND YSPEC_BUILD_EXAMPLES)
  file(GLOB_RECURSE YSPEC_EXAMPLE_FILES CONFIGURE_DEPENDS
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/examples
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.json
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.sh
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.pz
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/prem.*)
  set(YSPEC_EXAMPLE_COPIES)
  foreach(name IN LISTS YSPEC_EXAMPLE_FILES)
    set(file ${CMAKE_CURRENT_SOURCE_DIR}/examples/${name})
    set(copy ${CMAKE_CURRENT_BINARY_DIR}/examples/${name})
    add_custom_command(OUTPUT ${copy}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${file} ${copy}
      DEPENDS ${file}
      COMMENT "Copying examples/${name}")
    list(APPEND YSPEC_EXAMPLE_COPIES ${copy})
  endforeach()
  # And where the sources are, for the scripts to say where the Python
  # package, which draws the figures, is installed from.
  file(CONFIGURE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/examples/source.sh
    CONTENT "YSPEC_SOURCE_DIR=\"@PROJECT_SOURCE_DIR@\"\n" @ONLY)
  add_custom_target(examples ALL DEPENDS ${YSPEC_EXAMPLE_COPIES})
endif()

if(PROJECT_IS_TOP_LEVEL)
  option(YSPEC_BUILD_TESTS "Build the test suite" ON)
  if(YSPEC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
  endif()

  # The API documentation, on request: cmake --build build --target docs.
  find_package(Doxygen QUIET)
  if(Doxygen_FOUND)
    set(DOXYGEN_INPUT_DIR "${PROJECT_SOURCE_DIR}/YSpec")
    set(DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/docs")
    configure_file("${PROJECT_SOURCE_DIR}/docs/Doxyfile.in"
                   "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
    add_custom_target(docs
      COMMAND ${CMAKE_COMMAND} -E make_directory "${DOXYGEN_OUTPUT_DIR}"
      COMMAND ${DOXYGEN_EXECUTABLE} "${PROJECT_BINARY_DIR}/Doxyfile"
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
      COMMENT "Building the API documentation with Doxygen"
      VERBATIM)
  endif()
endif()
