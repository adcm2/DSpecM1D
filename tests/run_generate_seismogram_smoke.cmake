if(NOT DEFINED GENERATE_SEISMOGRAM_EXECUTABLE)
  message(FATAL_ERROR "GENERATE_SEISMOGRAM_EXECUTABLE is required")
endif()
if(NOT DEFINED GENERATE_SEISMOGRAM_PARAM_FILE)
  message(FATAL_ERROR "GENERATE_SEISMOGRAM_PARAM_FILE is required")
endif()
if(NOT DEFINED GENERATE_SEISMOGRAM_OUTPUT_FILE)
  message(FATAL_ERROR "GENERATE_SEISMOGRAM_OUTPUT_FILE is required")
endif()
if(NOT DEFINED GENERATE_SEISMOGRAM_WORKING_DIRECTORY)
  message(FATAL_ERROR "GENERATE_SEISMOGRAM_WORKING_DIRECTORY is required")
endif()

get_filename_component(_output_dir "${GENERATE_SEISMOGRAM_OUTPUT_FILE}" DIRECTORY)
file(MAKE_DIRECTORY "${_output_dir}")
file(REMOVE "${GENERATE_SEISMOGRAM_OUTPUT_FILE}")

execute_process(
  COMMAND "${GENERATE_SEISMOGRAM_EXECUTABLE}" "${GENERATE_SEISMOGRAM_PARAM_FILE}"
  WORKING_DIRECTORY "${GENERATE_SEISMOGRAM_WORKING_DIRECTORY}"
  RESULT_VARIABLE _generate_result
)

if(NOT _generate_result EQUAL 0)
  message(FATAL_ERROR
    "generate_seismogram smoke run failed with exit code ${_generate_result}")
endif()

if(NOT EXISTS "${GENERATE_SEISMOGRAM_OUTPUT_FILE}")
  message(FATAL_ERROR
    "generate_seismogram did not create expected output file:\n"
    "  ${GENERATE_SEISMOGRAM_OUTPUT_FILE}")
endif()

file(SIZE "${GENERATE_SEISMOGRAM_OUTPUT_FILE}" _output_size)
if(_output_size LESS 1)
  message(FATAL_ERROR
    "generate_seismogram created an empty output file:\n"
    "  ${GENERATE_SEISMOGRAM_OUTPUT_FILE}")
endif()
