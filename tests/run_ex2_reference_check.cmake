if(NOT DEFINED EX2_EXECUTABLE)
  message(FATAL_ERROR "EX2_EXECUTABLE is required")
endif()
if(NOT DEFINED EX2_TIME_EXPECTED)
  message(FATAL_ERROR "EX2_TIME_EXPECTED is required")
endif()
if(NOT DEFINED EX2_TIME_ACTUAL)
  message(FATAL_ERROR "EX2_TIME_ACTUAL is required")
endif()
if(NOT DEFINED EX2_FREQ_ACTUAL)
  message(FATAL_ERROR "EX2_FREQ_ACTUAL is required")
endif()

get_filename_component(_migration_dir "${EX2_TIME_ACTUAL}" DIRECTORY)
file(MAKE_DIRECTORY "${_migration_dir}")

set(ENV{DSPECM1D_EX2_TIME_OUT} "${EX2_TIME_ACTUAL}")
set(ENV{DSPECM1D_EX2_FREQ_OUT} "${EX2_FREQ_ACTUAL}")

execute_process(
  COMMAND "${EX2_EXECUTABLE}"
  RESULT_VARIABLE _ex2_result
)

if(NOT _ex2_result EQUAL 0)
  message(FATAL_ERROR "ex2 reference run failed with exit code ${_ex2_result}")
endif()

execute_process(
  COMMAND python3
          "${CMAKE_CURRENT_LIST_DIR}/compare_ex2_reference.py"
          "${EX2_TIME_EXPECTED}"
          "${EX2_TIME_ACTUAL}"
  RESULT_VARIABLE _compare_result
)

if(NOT _compare_result EQUAL 0)
  message(FATAL_ERROR
    "ex2 time-series output exceeded migration tolerance against the current checked-in reference:\n"
    "  expected: ${EX2_TIME_EXPECTED}\n"
    "  actual:   ${EX2_TIME_ACTUAL}")
endif()
