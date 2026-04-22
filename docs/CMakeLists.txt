set(DSPECM1D_DOCS_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs)
set(DSPECM1D_SITE_OUTPUT_DIR ${CMAKE_BINARY_DIR}/site)

find_package(Doxygen QUIET)
find_program(DSPECM1D_DOT_EXECUTABLE dot)

if(DSPECM1D_DOT_EXECUTABLE)
    set(DSPECM1D_DOXYGEN_HAVE_DOT YES)
else()
    set(DSPECM1D_DOXYGEN_HAVE_DOT NO)
endif()

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

find_package(Python3 COMPONENTS Interpreter QUIET)

if(Python3_Interpreter_FOUND)
    add_custom_target(preview-website
        COMMAND ${Python3_EXECUTABLE} -m http.server 8000 --directory ${DSPECM1D_SITE_OUTPUT_DIR}
        WORKING_DIRECTORY ${DSPECM1D_SITE_OUTPUT_DIR}
        COMMENT "Serving the DSpecM1D website locally at http://127.0.0.1:8000"
        VERBATIM
        USES_TERMINAL
    )
    add_dependencies(preview-website website)
endif()

if(DOXYGEN_FOUND)
    add_dependencies(website doxygen-docs)
    add_custom_command(TARGET website POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${DSPECM1D_DOCS_OUTPUT_DIR}/html
                ${DSPECM1D_SITE_OUTPUT_DIR}/api
        VERBATIM
    )
endif()
