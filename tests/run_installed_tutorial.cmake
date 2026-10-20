foreach(_required IN ITEMS CMAKE_COMMAND BUILD_DIR SOURCE_ROOT CMAKE_GENERATOR PYTHON_EXECUTABLE)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "${_required} is required")
  endif()
endforeach()
set(_prefix "${BUILD_DIR}/phase2-consumer-prefix")
set(_consumer_build "${BUILD_DIR}/phase2-consumer-build")
set(_input "${BUILD_DIR}/data/params/first_seismogram.txt")
set(_cli_output "${BUILD_DIR}/output/installed-check-cli.out")
set(_tutorial_output "${BUILD_DIR}/output/installed-check-consumer.out")
file(MAKE_DIRECTORY "${BUILD_DIR}/output")

execute_process(COMMAND "${CMAKE_COMMAND}" --install "${BUILD_DIR}" --prefix "${_prefix}"
  RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "Package installation failed: ${_result}")
endif()

set(_configure_args -S "${SOURCE_ROOT}/tutorials/first_seismogram" -B "${_consumer_build}"
  -DCMAKE_PREFIX_PATH=${_prefix} -DCMAKE_BUILD_TYPE=Release)
if(NOT "${CMAKE_GENERATOR}" STREQUAL "")
  list(APPEND _configure_args -G "${CMAKE_GENERATOR}")
endif()
foreach(_dep IN ITEMS
    "NUMERICCONCEPTS;${NUMERICCONCEPTS_SOURCE_DIR}"
    "FFTWPP;${FFTWPP_SOURCE_DIR}"
    "GAUSSQUAD;${GAUSSQUAD_SOURCE_DIR}"
    "GSHTRANS;${GSHTRANS_SOURCE_DIR}"
    "INTERPOLATION;${INTERPOLATION_SOURCE_DIR}")
  list(GET _dep 0 _name)
  list(GET _dep 1 _source)
  if(NOT "${_source}" STREQUAL "")
    list(APPEND _configure_args "-DFETCHCONTENT_SOURCE_DIR_${_name}=${_source}")
  endif()
endforeach()
execute_process(COMMAND "${CMAKE_COMMAND}" ${_configure_args} RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "Installed consumer configure failed: ${_result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" --build "${_consumer_build}" --parallel 2
  RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "Installed consumer build failed: ${_result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E env OMP_NUM_THREADS=2
  "${_consumer_build}/first_seismogram_consumer" "${_input}" "${_tutorial_output}"
  RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "Installed consumer run failed: ${_result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E env OMP_NUM_THREADS=2
  "${BUILD_DIR}/bin/generate_seismogram" "${_input}" RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "CLI comparison run failed: ${_result}")
endif()
set(_cli_output "${BUILD_DIR}/output/first_seismogram.out_t.out")
execute_process(COMMAND "${PYTHON_EXECUTABLE}"
  "${SOURCE_ROOT}/tests/compare_seismogram_outputs.py"
  "${_cli_output}" "${_tutorial_output}" RESULT_VARIABLE _compare_result)
if(NOT _compare_result EQUAL 0)
  message(FATAL_ERROR "Installed consumer output comparison failed: ${_compare_result}")
endif()
