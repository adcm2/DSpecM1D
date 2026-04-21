include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
  @ONLY
)

add_executable(dspecm1d_tests
  test_input_parser.cpp
  test_input_parameters_new.cpp
  test_output_writers.cpp
  test_reference_series_io.cpp
  test_signal_filtering.cpp
  test_spec_helpers.cpp
  test_sem_component.cpp
)

target_link_libraries(dspecm1d_tests PRIVATE
  DSpecM1D
  gtest_main
)

target_include_directories(dspecm1d_tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

include(GoogleTest)
gtest_discover_tests(
  dspecm1d_tests
  PROPERTIES
    LABELS "unit;component"
)
