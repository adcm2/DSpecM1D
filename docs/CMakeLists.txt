set(DSPECM1D_DOCS_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs)
set(DSPECM1D_SITE_OUTPUT_DIR ${CMAKE_BINARY_DIR}/site)

find_package(Doxygen QUIET)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    @ONLY
)

if(DOXYGEN_FOUND)
    add_custom_target(doxygen-docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Building DSpecM1D API documentation with Doxygen"
        VERBATIM
    )
endif()

add_custom_target(website
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${DSPECM1D_SITE_OUTPUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DSPECM1D_SITE_OUTPUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/site
            ${DSPECM1D_SITE_OUTPUT_DIR}
    COMMENT "Building the DSpecM1D static website"
    VERBATIM
)

if(DOXYGEN_FOUND)
    add_dependencies(website doxygen-docs)
    add_custom_command(TARGET website POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${DSPECM1D_DOCS_OUTPUT_DIR}/html
                ${DSPECM1D_SITE_OUTPUT_DIR}/api
        VERBATIM
    )
endif()
