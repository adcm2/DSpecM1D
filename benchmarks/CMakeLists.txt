# 1. Generate the config.h file inside the build/benchmarks directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in 
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

# 2. Grab sources and watch the folder for new additions
file(GLOB BENCHMARK_SOURCES CONFIGURE_DEPENDS "*.cpp")

# These retained solver-development programs include LAPACK/LAPACKE directly
# or require both backends.  Keep them out of an Eigen-only build so that
# configuration has no LAPACK dependency.
set(LAPACK_BENCHMARK_TARGETS
    eigen_pivot001_validation
    eigen_pivot0_validation
    eigen_pivot_error_characterization
    final_solver_benchmark
    lapack_band_concurrency_probe
    lapack_custom_backsolve_probe
    lapack_solve_allocation_probe
    lapacke_interface_probe
)

# 3. Loop through each file to create its own executable
foreach(SRC_FILE IN LISTS BENCHMARK_SOURCES)
    get_filename_component(TARGET_NAME ${SRC_FILE} NAME_WE)

    if(TARGET_NAME IN_LIST LAPACK_BENCHMARK_TARGETS
       AND NOT DSPECM1D_ENABLE_LAPACK_BAND_SOLVER)
        continue()
    endif()

    add_executable(${TARGET_NAME} ${SRC_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE DSpecM1D)
    target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

    if(TARGET CopyData)
        add_dependencies(${TARGET_NAME} CopyData)
    endif()

endforeach()

if(TARGET stage_profile_solver)
    target_compile_definitions(stage_profile_solver PRIVATE
        DSPECM1D_ENABLE_PROFILING)
endif()

if(TARGET ex1 AND TARGET ex2 AND TARGET ex3 AND TARGET ex4 AND TARGET ex5
   AND TARGET ex6 AND TARGET ex7)
    add_custom_target(paper_examples
        COMMENT "Build the paper-reproduction examples (ex1-ex7)"
        DEPENDS ex1 ex2 ex3 ex4 ex5 ex6 ex7
    )
endif()

if(DSPECM1D_ENABLE_SMOKE_TESTS AND TARGET ex1)
    add_test(
        NAME smoke_build_benchmark_ex1
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ex1
    )
    set_tests_properties(smoke_build_benchmark_ex1 PROPERTIES LABELS "dspecm1d_smoke")
endif()

if(DSPECM1D_ENABLE_MIGRATION_TESTS AND TARGET ex2)
    add_test(
        NAME migration_ex2_reference_check
        COMMAND ${CMAKE_COMMAND}
            -DEX2_EXECUTABLE=$<TARGET_FILE:ex2>
            -DEX2_TIME_EXPECTED=${CMAKE_SOURCE_DIR}/plotting/outputs/ex2_t.out
            -DEX2_TIME_ACTUAL=${CMAKE_BINARY_DIR}/migration_checks/ex2_t.out
            -DEX2_FREQ_ACTUAL=${CMAKE_BINARY_DIR}/migration_checks/ex2_w.out
            -P ${CMAKE_SOURCE_DIR}/tests/run_ex2_reference_check.cmake
    )
    set_tests_properties(migration_ex2_reference_check PROPERTIES
        LABELS "migration;reference"
    )
endif()

if(DSPECM1D_ENABLE_PAPER_VALIDATION AND TARGET paper_examples)
    add_test(
        NAME validate_build_paper_examples
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target paper_examples
    )
    set_tests_properties(validate_build_paper_examples PROPERTIES
        LABELS "validation;paper"
    )
endif()
