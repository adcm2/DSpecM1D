file(GLOB TOOL_SOURCES CONFIGURE_DEPENDS "*.cpp")

foreach(SRC_FILE IN LISTS TOOL_SOURCES)
    get_filename_component(TARGET_NAME ${SRC_FILE} NAME_WE)
    add_executable(${TARGET_NAME} ${SRC_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE DSpecM1D)

    if(TARGET CopyData)
        add_dependencies(${TARGET_NAME} CopyData)
    endif()
endforeach()

if(DSPECM1D_ENABLE_SMOKE_TESTS AND TARGET generate_seismogram)
    add_test(
        NAME smoke_generate_seismogram
        COMMAND ${CMAKE_COMMAND}
            -DGENERATE_SEISMOGRAM_EXECUTABLE=$<TARGET_FILE:generate_seismogram>
            -DGENERATE_SEISMOGRAM_PARAM_FILE=${CMAKE_BINARY_DIR}/data/params/generate_seismogram_smoke.txt
            -DGENERATE_SEISMOGRAM_OUTPUT_FILE=${CMAKE_BINARY_DIR}/output/generate_seismogram_smoke_t.out
            -DGENERATE_SEISMOGRAM_WORKING_DIRECTORY=${CMAKE_BINARY_DIR}
            -P ${CMAKE_SOURCE_DIR}/tests/run_generate_seismogram_smoke.cmake
    )
    set_tests_properties(smoke_generate_seismogram PROPERTIES
        LABELS "dspecm1d_smoke"
        RUN_SERIAL TRUE
        TIMEOUT 600
        ENVIRONMENT "OMP_NUM_THREADS=2"
    )
endif()
