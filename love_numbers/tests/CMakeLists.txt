add_executable(dspecm1d_love_numbers_tests
    test_love_numbers.cpp
)

add_executable(dspecm1d_love_numbers_radial_state_tests
    test_radial_state.cpp
)

add_executable(dspecm1d_love_numbers_dof_map_tests
    test_love_dof_map.cpp
)

add_executable(dspecm1d_love_numbers_static_operator_tests
    test_static_operator.cpp
)

add_executable(dspecm1d_love_numbers_static_fluid_tests
    test_static_fluids.cpp
)

add_executable(dspecm1d_love_numbers_static_solve_tests
    test_static_solve.cpp
)

add_executable(dspecm1d_love_numbers_static_radial_tests
    test_static_radial.cpp
)

add_executable(dspecm1d_love_numbers_static_degree_one_tests
    test_static_degree_one.cpp
)

add_executable(dspecm1d_love_numbers_public_private_tests
    test_public_private.cpp
)

add_executable(dspecm1d_love_numbers_cli_tests
    test_cli.cpp
)

if(DSPECM1D_ENABLE_PAPER_VALIDATION)
    add_executable(dspecm1d_love_numbers_reference_validation
        test_reference_validation.cpp
    )
endif()

target_link_libraries(dspecm1d_love_numbers_tests
    PRIVATE
        DSpecM1D::LoveNumbers
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_radial_state_tests
    PRIVATE
        DSpecM1D
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_dof_map_tests
    PRIVATE
        DSpecM1D
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_static_operator_tests
    PRIVATE
        DSpecM1D
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_static_fluid_tests
    PRIVATE
        DSpecM1D
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_static_solve_tests
    PRIVATE
        DSpecM1D
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_static_radial_tests
    PRIVATE
        DSpecM1D
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_static_degree_one_tests
    PRIVATE
        DSpecM1D
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_public_private_tests
    PRIVATE
        DSpecM1D
        gtest_main
)

target_link_libraries(dspecm1d_love_numbers_cli_tests
    PRIVATE
        DSpecM1D::LoveNumbers
        gtest_main
)

add_dependencies(dspecm1d_love_numbers_cli_tests dspecm1d-love)

if(DSPECM1D_ENABLE_PAPER_VALIDATION)
    target_link_libraries(dspecm1d_love_numbers_reference_validation
        PRIVATE
            DSpecM1D::LoveNumbers
            gtest_main
    )
endif()

target_compile_definitions(dspecm1d_love_numbers_radial_state_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.noatten.txt"
        DSPECM1D_LOVE_NUMBERS_OCEAN_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.noatten.txt"
        DSPECM1D_LOVE_NUMBERS_HOMOGENEOUS_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/../validation/data/homogeneous_sphere.dspec"
        DSPECM1D_LOVE_NUMBERS_CONTROLLED_SFS_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/../validation/data/solid_fluid_solid.dspec"
        DSPECM1D_LOVE_NUMBERS_CONTROLLED_CENTRAL_FLUID_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/../validation/data/central_fluid.dspec"
        DSPECM1D_LOVE_NUMBERS_SFS_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/internal_fluid_sfs.txt"
)

target_compile_definitions(dspecm1d_love_numbers_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_isotropic.txt"
        DSPECM1D_LOVE_NUMBERS_TI_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_ti.txt"
        DSPECM1D_LOVE_NUMBERS_SFS_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/internal_fluid_sfs.txt"
        DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.noatten.txt"
        DSPECM1D_LOVE_NUMBERS_OCEAN_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.noatten.txt"
)

target_compile_definitions(dspecm1d_love_numbers_dof_map_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.noatten.txt"
)

target_compile_definitions(dspecm1d_love_numbers_static_operator_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_isotropic.txt"
        DSPECM1D_LOVE_NUMBERS_TI_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_ti.txt"
        DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.noatten.txt"
)

target_compile_definitions(dspecm1d_love_numbers_static_fluid_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_SFS_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/internal_fluid_sfs.txt"
        DSPECM1D_LOVE_NUMBERS_CENTRAL_FLUID_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/central_constant_fluid.txt"
        DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.noatten.txt"
)

target_compile_definitions(dspecm1d_love_numbers_static_solve_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_isotropic.txt"
        DSPECM1D_LOVE_NUMBERS_TI_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_ti.txt"
        DSPECM1D_LOVE_NUMBERS_SFS_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/internal_fluid_sfs.txt"
        DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.noatten.txt"
)

target_compile_definitions(dspecm1d_love_numbers_static_radial_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_isotropic.txt"
        DSPECM1D_LOVE_NUMBERS_TI_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_ti.txt"
        DSPECM1D_LOVE_NUMBERS_SFS_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/internal_fluid_sfs.txt"
        DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.noatten.txt"
)

target_compile_definitions(dspecm1d_love_numbers_static_degree_one_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_isotropic.txt"
        DSPECM1D_LOVE_NUMBERS_TI_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_ti.txt"
        DSPECM1D_LOVE_NUMBERS_SFS_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/internal_fluid_sfs.txt"
        DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.noatten.txt"
)

target_compile_definitions(dspecm1d_love_numbers_public_private_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_TI_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_ti.txt"
)

target_compile_definitions(dspecm1d_love_numbers_cli_tests
    PRIVATE
        DSPECM1D_LOVE_NUMBERS_CLI="$<TARGET_FILE:dspecm1d-love>"
        DSPECM1D_LOVE_NUMBERS_CLI_MODEL="${CMAKE_CURRENT_SOURCE_DIR}/data/all_solid_isotropic.txt"
        DSPECM1D_LOVE_NUMBERS_CLI_OCEAN_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.noatten.txt"
)

if(DSPECM1D_ENABLE_PAPER_VALIDATION)
    target_compile_definitions(dspecm1d_love_numbers_reference_validation
        PRIVATE
            DSPECM1D_LOVE_NUMBERS_REFERENCE="${CMAKE_CURRENT_SOURCE_DIR}/data/prem_love_reference.dat"
            DSPECM1D_LOVE_NUMBERS_REFERENCE_MODEL="${PROJECT_SOURCE_DIR}/data/models/prem.200.no.txt"
    )
endif()

target_include_directories(dspecm1d_love_numbers_radial_state_tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

target_include_directories(dspecm1d_love_numbers_dof_map_tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

target_include_directories(dspecm1d_love_numbers_static_operator_tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

target_include_directories(dspecm1d_love_numbers_static_fluid_tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

target_include_directories(dspecm1d_love_numbers_static_solve_tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

target_include_directories(dspecm1d_love_numbers_static_radial_tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

target_include_directories(dspecm1d_love_numbers_static_degree_one_tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

target_include_directories(dspecm1d_love_numbers_public_private_tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
        ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

include(GoogleTest)
gtest_discover_tests(
    dspecm1d_love_numbers_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_radial_state_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_dof_map_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_static_operator_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_static_fluid_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_static_solve_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_static_radial_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_static_degree_one_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_public_private_tests
    PROPERTIES
        LABELS "unit;component"
)

gtest_discover_tests(
    dspecm1d_love_numbers_cli_tests
    PROPERTIES
        LABELS "unit;component"
)

if(DSPECM1D_ENABLE_PAPER_VALIDATION)
    gtest_discover_tests(
        dspecm1d_love_numbers_reference_validation
        PROPERTIES
            LABELS "paper-validation"
    )
endif()
