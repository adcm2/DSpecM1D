add_library(DSpecM1DLoveNumbers
    src/LoveNumbers.cpp
    src/LoveDofMap.h
    src/StaticForcing.h
    src/StaticOperator.h
    src/StaticSolve.h
)
add_library(DSpecM1D::LoveNumbers ALIAS DSpecM1DLoveNumbers)

target_link_libraries(DSpecM1DLoveNumbers
    PUBLIC DSpecM1D
)

add_executable(dspecm1d-love
    app/dspecm1d-love.cpp
)
target_link_libraries(dspecm1d-love
    PRIVATE DSpecM1D::LoveNumbers
)

target_include_directories(DSpecM1DLoveNumbers
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

if(PROJECT_IS_TOP_LEVEL AND DSPECM1D_BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(PROJECT_IS_TOP_LEVEL AND DSPECM1D_ENABLE_PAPER_VALIDATION)
    enable_language(Fortran)
    add_subdirectory(validation)
endif()
