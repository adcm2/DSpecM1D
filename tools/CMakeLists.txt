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
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_test(
        NAME smoke_generate_seismogram
        COMMAND ${CMAKE_COMMAND}
            -DGENERATE_SEISMOGRAM_EXECUTABLE=$<TARGET_FILE:generate_seismogram>
            -DGENERATE_SEISMOGRAM_PARAM_FILE=${CMAKE_BINARY_DIR}/data/params/first_seismogram.txt
            -DGENERATE_SEISMOGRAM_OUTPUT_FILE=${CMAKE_BINARY_DIR}/output/first_seismogram.out_t.out
            -DGENERATE_SEISMOGRAM_WORKING_DIRECTORY=${CMAKE_BINARY_DIR}
            -DGENERATE_SEISMOGRAM_CHECKER=${CMAKE_SOURCE_DIR}/tests/check_first_seismogram.py
            -DGENERATE_SEISMOGRAM_PYTHON=${Python3_EXECUTABLE}
            -P ${CMAKE_SOURCE_DIR}/tests/run_generate_seismogram_smoke.cmake
    )
    set_tests_properties(smoke_generate_seismogram PROPERTIES
        LABELS "dspecm1d_smoke"
        RUN_SERIAL TRUE
        TIMEOUT 600
        ENVIRONMENT "OMP_NUM_THREADS=2"
    )
endif()

if(DSPECM1D_ENABLE_SMOKE_TESTS AND TARGET generate_seismogram AND TARGET first_seismogram)
    FetchContent_GetProperties(NumericConcepts SOURCE_DIR _dspecm1d_numericconcepts_source)
    FetchContent_GetProperties(FFTWpp SOURCE_DIR _dspecm1d_fftwpp_source)
    FetchContent_GetProperties(GaussQuad SOURCE_DIR _dspecm1d_gaussquad_source)
    FetchContent_GetProperties(GSHTRANS SOURCE_DIR _dspecm1d_gshtrans_source)
    FetchContent_GetProperties(Interpolation SOURCE_DIR _dspecm1d_interpolation_source)
    add_test(
        NAME smoke_installed_first_seismogram
        COMMAND ${CMAKE_COMMAND}
            -DCMAKE_COMMAND=${CMAKE_COMMAND}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/tutorials/first_seismogram
            -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -DSOURCE_ROOT=${CMAKE_SOURCE_DIR}
            -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
            -DPYTHON_EXECUTABLE=${Python3_EXECUTABLE}
            -DNUMERICCONCEPTS_SOURCE_DIR=${_dspecm1d_numericconcepts_source}
            -DFFTWPP_SOURCE_DIR=${_dspecm1d_fftwpp_source}
            -DGAUSSQUAD_SOURCE_DIR=${_dspecm1d_gaussquad_source}
            -DGSHTRANS_SOURCE_DIR=${_dspecm1d_gshtrans_source}
            -DINTERPOLATION_SOURCE_DIR=${_dspecm1d_interpolation_source}
            -P ${CMAKE_SOURCE_DIR}/tests/run_installed_tutorial.cmake
    )
    set_tests_properties(smoke_installed_first_seismogram PROPERTIES
        LABELS "dspecm1d_smoke;installed-package;tutorial"
        RUN_SERIAL TRUE
        TIMEOUT 1200
        ENVIRONMENT "OMP_NUM_THREADS=2"
    )
endif()
