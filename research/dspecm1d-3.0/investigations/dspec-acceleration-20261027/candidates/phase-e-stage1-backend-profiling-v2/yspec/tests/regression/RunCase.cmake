# Runs one regression case: copies its input and model into a fresh working
# directory, runs yspec there, and compares each output with the reference of
# the same name. Invoked by ctest with -P; see tests/CMakeLists.txt.

file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR})
file(COPY ${CASE_DIR}/input.json DESTINATION ${WORK_DIR})
file(READ ${CASE_DIR}/input.json input)
string(JSON model GET "${input}" model)
file(COPY ${EXAMPLES_DIR}/${model} DESTINATION ${WORK_DIR})

execute_process(COMMAND ${YSPEC} input.json
                WORKING_DIRECTORY ${WORK_DIR}
                ERROR_FILE ${WORK_DIR}/yspec.log
                RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "yspec exited with status ${status}; see ${WORK_DIR}/yspec.log")
endif()

file(GLOB references ${CASE_DIR}/reference/*)
if(NOT references)
  message(FATAL_ERROR "no reference files in ${CASE_DIR}/reference")
endif()
set(failed FALSE)
foreach(reference ${references})
  get_filename_component(name ${reference} NAME)
  execute_process(COMMAND ${COMPARE} ${WORK_DIR}/${name} ${reference} ${TOLERANCE}
                  RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    set(failed TRUE)
  endif()
endforeach()
if(failed)
  message(FATAL_ERROR "outputs differ from the references")
endif()
