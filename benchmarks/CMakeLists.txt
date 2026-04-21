# 1. Generate the config.h file inside the build/benchmarks directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in 
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

# 2. Grab sources and watch the folder for new additions
file(GLOB BENCHMARK_SOURCES CONFIGURE_DEPENDS "*.cpp")

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

if(DSPECM1D_ENABLE_SMOKE_TESTS AND TARGET ex1)
    add_test(
        NAME smoke_build_benchmark_ex1
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ex1
    )
    set_tests_properties(smoke_build_benchmark_ex1 PROPERTIES LABELS "smoke")
endif()