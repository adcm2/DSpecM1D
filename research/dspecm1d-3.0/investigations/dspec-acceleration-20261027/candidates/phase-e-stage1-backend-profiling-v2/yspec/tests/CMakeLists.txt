FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.17.0
  FIND_PACKAGE_ARGS NAMES GTest)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

# Tests read the models in examples/ and the reference outputs in tests/data/.
set(YSPEC_TEST_DEFINITIONS
    YSPEC_EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/examples"
    YSPEC_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

# Tests of the internals, in src/, link YSpecInternal as well as the library.
function(yspec_add_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE YSpec::YSpec YSpecInternal
                                        GTest::gtest_main)
  target_compile_definitions(${name} PRIVATE ${YSPEC_TEST_DEFINITIONS})
  target_compile_options(${name} PRIVATE ${YSPEC_WARNINGS})
  gtest_discover_tests(${name} DISCOVERY_MODE PRE_TEST)
endfunction()

# Tests of a full workflow, labelled slow: everyday runs can leave them out
# with ctest -LE slow, as the Debug builds of CI do. Each has a quick
# counterpart that tests the pieces.
function(yspec_add_slow_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE YSpec::YSpec YSpecInternal
                                        GTest::gtest_main)
  target_compile_definitions(${name} PRIVATE ${YSPEC_TEST_DEFINITIONS})
  target_compile_options(${name} PRIVATE ${YSPEC_WARNINGS})
  gtest_discover_tests(${name} DISCOVERY_MODE PRE_TEST PROPERTIES LABELS slow)
endfunction()

yspec_add_test(TestModel)
yspec_add_test(TestNumerics)
yspec_add_test(TestSolver)
yspec_add_test(TestBoundaries)
yspec_add_test(TestReceivers)
yspec_add_test(TestInput)
yspec_add_test(TestInstrument)
yspec_add_test(TestNormalModes)
yspec_add_test(TestHomogeneousSphere)
yspec_add_test(TestSynthetics)
yspec_add_test(TestStrain)
yspec_add_test(TestKernels)
yspec_add_test(TestKernelParts)
yspec_add_test(TestVolumeKernels)
yspec_add_test(TestSnapshots)
yspec_add_test(TestRun)
yspec_add_slow_test(TestKernelsWorkflow)
yspec_add_slow_test(TestStrainWorkflow)

# End-to-end cases, compared with seismograms from the Fortran YSpec (made by
# tests/fortran/make_regression_reference.sh). Differences are relative to the
# peak of each component.
add_executable(CompareSeismograms regression/CompareSeismograms.cpp)
set(YSPEC_REGRESSION_TOLERANCE 5e-4)

function(yspec_add_case name)
  add_test(NAME regression.${name}
           COMMAND ${CMAKE_COMMAND}
                   -DYSPEC=$<TARGET_FILE:yspec>
                   -DCOMPARE=$<TARGET_FILE:CompareSeismograms>
                   -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/regression/${name}
                   -DEXAMPLES_DIR=${PROJECT_SOURCE_DIR}/examples
                   -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/regression/${name}
                   -DTOLERANCE=${YSPEC_REGRESSION_TOLERANCE}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/regression/RunCase.cmake)
endfunction()

# Full gravity, solid inner core, attenuation, displacement; receivers from 5
# to 175 degrees and at the north pole.
yspec_add_case(selfgrav_disp)
# Cowling, attenuation, acceleration with potential and tilt corrections; a
# source in the low-Q zone and a receiver at the south pole.
yspec_add_case(cowling_accel)
# No gravity, fluid inner core, no attenuation, velocity; receiver below the
# source.
yspec_add_case(fic_nograv_vel)
# Radial modes only, from an isotropic source; acceleration with the free-air
# correction, which is large near 0S0.
yspec_add_case(radial_accel)

# The example in README.md.
add_executable(library_example ${PROJECT_SOURCE_DIR}/examples/library.cpp)
target_link_libraries(library_example PRIVATE YSpec::YSpec)
target_compile_options(library_example PRIVATE ${YSPEC_WARNINGS})
target_compile_definitions(library_example PRIVATE ${YSPEC_TEST_DEFINITIONS})
add_test(NAME example.library COMMAND library_example)
