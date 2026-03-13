cmake_minimum_required(VERSION 3.24)

# Ensure data is copied before the test runs
add_test(
    NAME    setup_data
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target CopyData
)

add_test(
    NAME    smoke_ex3
    COMMAND $<TARGET_FILE:ex3>
            ${CMAKE_BINARY_DIR}/data/params/ex3.txt
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/data
)

set_tests_properties(smoke_ex3 PROPERTIES
    TIMEOUT  300
    DEPENDS  setup_data
)
