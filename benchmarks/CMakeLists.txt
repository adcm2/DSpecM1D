# add_executable(ex1 ex1.cpp)
# target_link_libraries(ex1 PRIVATE DSpecM1D)

# add_executable(ex2 ex2.cpp)
# target_link_libraries(ex2 PRIVATE DSpecM1D)

# add_executable(ex3 ex3.cpp)
# target_link_libraries(ex3 PRIVATE DSpecM1D)

# add_executable(ex4 ex4.cpp)
# target_link_libraries(ex4 PRIVATE DSpecM1D)

# add_executable(ex5 ex5.cpp)
# target_link_libraries(ex5 PRIVATE DSpecM1D)

# add_executable(ex6 ex6.cpp)
# target_link_libraries(ex6 PRIVATE DSpecM1D)

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