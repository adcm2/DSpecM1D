include(FetchContent)

set(DSPECM1D_GIA3D_SOURCE_DIR "" CACHE PATH
    "Optional local da380/gia3D checkout for Love-number validation"
)
set(DSPECM1D_GIA3D_CORE_SOURCE_DIR "" CACHE PATH
    "Optional local da380/core checkout for Love-number validation"
)

set(gia3d_revision e6a8230fc792416087e530726e3f9f9b41ee0506)
set(gia3d_core_revision eb1d36e418ab0ac9e0b34abb717e6d9b2e863ee5)

if(DSPECM1D_GIA3D_SOURCE_DIR)
    set(gia3d_source_dir "${DSPECM1D_GIA3D_SOURCE_DIR}")
else()
    FetchContent_Declare(
        gia3d_validation_source
        GIT_REPOSITORY https://github.com/da380/gia3D.git
        GIT_TAG ${gia3d_revision}
        GIT_SHALLOW FALSE
    )
    FetchContent_MakeAvailable(gia3d_validation_source)
    set(gia3d_source_dir "${gia3d_validation_source_SOURCE_DIR}")
endif()

if(DSPECM1D_GIA3D_CORE_SOURCE_DIR)
    set(gia3d_core_source_dir "${DSPECM1D_GIA3D_CORE_SOURCE_DIR}")
else()
    FetchContent_Declare(
        gia3d_core_validation_source
        GIT_REPOSITORY https://github.com/da380/core.git
        GIT_TAG ${gia3d_core_revision}
        GIT_SHALLOW FALSE
    )
    FetchContent_MakeAvailable(gia3d_core_validation_source)
    set(gia3d_core_source_dir "${gia3d_core_validation_source_SOURCE_DIR}")
endif()

find_package(Git REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_program(PATCH_EXECUTABLE patch REQUIRED)

function(verify_validation_revision source_dir expected_revision label)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} -C "${source_dir}" rev-parse HEAD
        RESULT_VARIABLE revision_result
        OUTPUT_VARIABLE actual_revision
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT revision_result EQUAL 0 OR
       NOT actual_revision STREQUAL expected_revision)
        message(FATAL_ERROR
            "${label} must be checkout ${expected_revision}; got "
            "${actual_revision} from ${source_dir}"
        )
    endif()
endfunction()

verify_validation_revision(
    "${gia3d_source_dir}" "${gia3d_revision}" "gia3D validation source"
)
verify_validation_revision(
    "${gia3d_core_source_dir}" "${gia3d_core_revision}"
    "gia3D core validation source"
)

set(patched_source_dir "${CMAKE_CURRENT_BINARY_DIR}/patched-gia3d")
file(MAKE_DIRECTORY "${patched_source_dir}/src")
configure_file(
    "${gia3d_source_dir}/src/module_matrix.f90"
    "${patched_source_dir}/src/module_matrix.f90"
    COPYONLY
)
configure_file(
    "${gia3d_source_dir}/src/module_force.f90"
    "${patched_source_dir}/src/module_force.f90"
    COPYONLY
)
execute_process(
    COMMAND
        ${PATCH_EXECUTABLE} --batch --forward -p1
        --input=${CMAKE_CURRENT_SOURCE_DIR}/gia3d/patches/initialize_fluid_radius.patch
    WORKING_DIRECTORY "${patched_source_dir}"
    RESULT_VARIABLE patch_result
    OUTPUT_VARIABLE patch_output
    ERROR_VARIABLE patch_error
)
if(NOT patch_result EQUAL 0)
    message(FATAL_ERROR
        "Could not patch gia3D module_matrix.f90:\n"
        "${patch_output}${patch_error}"
    )
endif()
execute_process(
    COMMAND
        ${PATCH_EXECUTABLE} --batch --forward -p1
        --input=${CMAKE_CURRENT_SOURCE_DIR}/gia3d/patches/safe_central_fluid_tide.patch
    WORKING_DIRECTORY "${patched_source_dir}"
    RESULT_VARIABLE force_patch_result
    OUTPUT_VARIABLE force_patch_output
    ERROR_VARIABLE force_patch_error
)
if(NOT force_patch_result EQUAL 0)
    message(FATAL_ERROR
        "Could not patch gia3D module_force.f90:\n"
        "${force_patch_output}${force_patch_error}"
    )
endif()

set(gia3d_core_sources
    "${gia3d_core_source_dir}/src/module_constants.f90"
    "${gia3d_core_source_dir}/src/module_physical_constants.f90"
    "${gia3d_core_source_dir}/src/module_error.f90"
    "${gia3d_core_source_dir}/src/module_util.f90"
    "${gia3d_core_source_dir}/src/module_interp.f90"
    "${gia3d_core_source_dir}/src/module_special_functions.f90"
    "${gia3d_core_source_dir}/src/module_quadrature.f90"
    "${gia3d_core_source_dir}/src/module_spherical_model.f90"
    "${gia3d_core_source_dir}/src/module_DECK.f90"
    "${gia3d_core_source_dir}/src/module_PREM.f90"
    "${gia3d_core_source_dir}/src/module_spherical_mesh.f90"
)

add_library(dspecm1d_gia3d_reference_modules STATIC
    ${gia3d_core_sources}
    "${patched_source_dir}/src/module_matrix.f90"
    "${patched_source_dir}/src/module_force.f90"
    gia3d/full_domain_matrix.f90
)
set(patched_module_dir "${CMAKE_CURRENT_BINARY_DIR}/modules")
set_target_properties(dspecm1d_gia3d_reference_modules PROPERTIES
    Fortran_MODULE_DIRECTORY "${patched_module_dir}"
)
target_include_directories(dspecm1d_gia3d_reference_modules
    PUBLIC "${patched_module_dir}"
)

add_executable(dspecm1d_gia3d_direct_reference
    gia3d/direct_love_numbers.f90
)
target_link_libraries(dspecm1d_gia3d_direct_reference PRIVATE
    dspecm1d_gia3d_reference_modules
    LAPACK::LAPACK
    BLAS::BLAS
)

add_executable(dspecm1d_gia3d_controlled_reference
    gia3d/controlled_love_numbers.f90
)
target_link_libraries(dspecm1d_gia3d_controlled_reference PRIVATE
    dspecm1d_gia3d_reference_modules
    LAPACK::LAPACK
    BLAS::BLAS
)

add_executable(dspecm1d_controlled_love_numbers
    controlled_love_numbers.cpp
)
target_include_directories(dspecm1d_controlled_love_numbers PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../include"
)
target_link_libraries(dspecm1d_controlled_love_numbers PRIVATE
    DSpecM1D
)

add_library(dspecm1d_gia3d_unpatched_modules STATIC EXCLUDE_FROM_ALL
    ${gia3d_core_sources}
    "${gia3d_source_dir}/src/module_matrix.f90"
    "${gia3d_source_dir}/src/module_force.f90"
)
set(unpatched_module_dir "${CMAKE_CURRENT_BINARY_DIR}/unpatched-modules")
set_target_properties(dspecm1d_gia3d_unpatched_modules PROPERTIES
    Fortran_MODULE_DIRECTORY "${unpatched_module_dir}"
)
target_include_directories(dspecm1d_gia3d_unpatched_modules
    PUBLIC "${unpatched_module_dir}"
)

add_executable(
    dspecm1d_gia3d_direct_reference_unpatched_checked EXCLUDE_FROM_ALL
    gia3d/direct_love_numbers.f90
)
target_link_libraries(
    dspecm1d_gia3d_direct_reference_unpatched_checked PRIVATE
    dspecm1d_gia3d_unpatched_modules
    LAPACK::LAPACK
    BLAS::BLAS
)

set(runtime_check_flags
    -O0
    -g
    -fcheck=all
    -finit-real=snan
    -fbacktrace
)
target_compile_options(dspecm1d_gia3d_unpatched_modules PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:${runtime_check_flags}>
)
target_compile_options(
    dspecm1d_gia3d_direct_reference_unpatched_checked PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:${runtime_check_flags}>
)

string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_upper)
string(REPLACE ";" " " blas_libraries "${BLAS_LIBRARIES}")
string(REPLACE ";" " " lapack_libraries "${LAPACK_LIBRARIES}")
file(GENERATE
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/gia3d-build-provenance.txt"
    CONTENT
"gia3D revision: ${gia3d_revision}
core revision: ${gia3d_core_revision}
gia3D source: ${gia3d_source_dir}
core source: ${gia3d_core_source_dir}
Fortran compiler: ${CMAKE_Fortran_COMPILER}
Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}
Fortran compiler version: ${CMAKE_Fortran_COMPILER_VERSION}
Build type: ${CMAKE_BUILD_TYPE}
Fortran flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${build_type_upper}}
Unpatched runtime-check flags: ${runtime_check_flags}
BLAS libraries: ${blas_libraries}
LAPACK libraries: ${lapack_libraries}
Matrix patch: ${CMAKE_CURRENT_SOURCE_DIR}/gia3d/patches/initialize_fluid_radius.patch
Tide patch: ${CMAKE_CURRENT_SOURCE_DIR}/gia3d/patches/safe_central_fluid_tide.patch
"
)

add_test(
    NAME LoveNumbersGia3DReference.ReportsDirectComparison
    COMMAND
        ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/gia3d/compare_reference.py"
        "$<TARGET_FILE:dspecm1d_gia3d_direct_reference>"
        "${CMAKE_CURRENT_SOURCE_DIR}/../tests/data/prem_love_reference.dat"
        1 2 10
)
set_tests_properties(
    LoveNumbersGia3DReference.ReportsDirectComparison
    PROPERTIES LABELS "paper-validation"
)

add_test(
    NAME LoveNumbersGia3DReference.ReportsControlledSphereComparison
    COMMAND
        ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/compare_controlled_sphere.py"
        "$<TARGET_FILE:dspecm1d_controlled_love_numbers>"
        "$<TARGET_FILE:dspecm1d_gia3d_controlled_reference>"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/homogeneous_sphere.dspec"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/homogeneous_sphere.gia3d"
)
set_tests_properties(
    LoveNumbersGia3DReference.ReportsControlledSphereComparison
    PROPERTIES LABELS "paper-validation"
)

add_test(
    NAME LoveNumbersGia3DReference.ReportsSolidFluidSolidComparison
    COMMAND
        ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/compare_controlled_sphere.py"
        "$<TARGET_FILE:dspecm1d_controlled_love_numbers>"
        "$<TARGET_FILE:dspecm1d_gia3d_controlled_reference>"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/solid_fluid_solid.dspec"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/solid_fluid_solid.gia3d"
        solid-fluid-solid
)
set_tests_properties(
    LoveNumbersGia3DReference.ReportsSolidFluidSolidComparison
    PROPERTIES LABELS "paper-validation"
)

add_test(
    NAME LoveNumbersGia3DReference.ReportsCentralFluidComparison
    COMMAND
        ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/compare_controlled_sphere.py"
        "$<TARGET_FILE:dspecm1d_controlled_love_numbers>"
        "$<TARGET_FILE:dspecm1d_gia3d_controlled_reference>"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/central_fluid.dspec"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/central_fluid.gia3d"
        central-fluid
)
set_tests_properties(
    LoveNumbersGia3DReference.ReportsCentralFluidComparison
    PROPERTIES LABELS "paper-validation"
)
