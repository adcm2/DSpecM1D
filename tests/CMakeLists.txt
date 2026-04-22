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

include(GoogleTest)

function(dspecm1d_add_test_exe target)
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} PRIVATE DSpecM1D gtest_main)
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  gtest_discover_tests(
    ${target}
    PROPERTIES
      LABELS "unit;component"
  )
endfunction()

dspecm1d_add_test_exe(dspecm1d_parser_tests
  test_input_parser.cpp
  test_reference_series_io.cpp
  test_spec_helpers.cpp
)

dspecm1d_add_test_exe(dspecm1d_input_context_tests
  test_input_parameters_new.cpp
)

dspecm1d_add_test_exe(dspecm1d_output_writer_tests
  test_output_writers.cpp
)

dspecm1d_add_test_exe(dspecm1d_signal_filter_tests
  test_signal_filtering.cpp
)

dspecm1d_add_test_exe(dspecm1d_sem_component_tests
  test_sem_component.cpp
)
