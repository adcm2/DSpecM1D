# 1. Generate the config.h file inside the build/benchmarks directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in 
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

file(GLOB BENCHMARK_SOURCES "*.cpp")

# 3. Loop through each file to create its own executable
foreach(SRC_FILE IN LISTS BENCHMARK_SOURCES)
    get_filename_component(TARGET_NAME ${SRC_FILE} NAME_WE)
    add_executable(${TARGET_NAME} ${SRC_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE DSpecM1D)
    target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

    
    if(TARGET CopyData)
        add_dependencies(${TARGET_NAME} CopyData)
    endif()

endforeach()

if(DSPECM1D_ENABLE_SMOKE_TESTS AND TARGET t1)
    add_test(NAME smoke_tutorial_t1 COMMAND $<TARGET_FILE:t1>)
    set_tests_properties(smoke_tutorial_t1 PROPERTIES
        LABELS "dspecm1d_smoke;tutorial"
        RUN_SERIAL TRUE
        TIMEOUT 1800
        ENVIRONMENT "OMP_NUM_THREADS=2"
    )
endif()

if(DSPECM1D_ENABLE_SMOKE_TESTS AND TARGET first_seismogram)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_test(
        NAME smoke_tutorial_first_seismogram
        COMMAND ${CMAKE_COMMAND}
            -DGENERATE_SEISMOGRAM_EXECUTABLE=$<TARGET_FILE:first_seismogram>
            -DGENERATE_SEISMOGRAM_PARAM_FILE=${CMAKE_BINARY_DIR}/data/params/first_seismogram.txt
            -DGENERATE_SEISMOGRAM_OUTPUT_FILE=${CMAKE_BINARY_DIR}/output/tutorial_first_seismogram_t.out
            -DGENERATE_SEISMOGRAM_WORKING_DIRECTORY=${CMAKE_BINARY_DIR}
            -DGENERATE_SEISMOGRAM_CHECKER=${CMAKE_SOURCE_DIR}/tests/check_first_seismogram.py
            -DGENERATE_SEISMOGRAM_PYTHON=${Python3_EXECUTABLE}
            -P ${CMAKE_SOURCE_DIR}/tests/run_tutorial_seismogram_smoke.cmake
    )
    set_tests_properties(smoke_tutorial_first_seismogram PROPERTIES
        LABELS "dspecm1d_smoke;tutorial"
        RUN_SERIAL TRUE
        TIMEOUT 600
        ENVIRONMENT "OMP_NUM_THREADS=2"
    )
endif()
