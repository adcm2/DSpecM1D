include(FetchContent)

set(DSPECM1D_GIA3D_SOURCE_DIR "" CACHE PATH
    "Optional local da380/gia3D checkout for Love-number validation"
)
set(DSPECM1D_GIA3D_CORE_SOURCE_DIR "" CACHE PATH
    "Optional local da380/core checkout for Love-number validation"
)
set(DSPECM1D_ELLN_ARCHIVE "" CACHE FILEPATH
    "Optional official ggx444_supp.zip archive for ELLN validation"
)
set(DSPECM1D_ELLN_SOURCE_DIR "" CACHE PATH
    "Optional extracted official ELLN source directory"
)

if(DSPECM1D_ELLN_ARCHIVE AND DSPECM1D_ELLN_SOURCE_DIR)
    message(FATAL_ERROR
        "Set only one of DSPECM1D_ELLN_ARCHIVE and "
        "DSPECM1D_ELLN_SOURCE_DIR"
    )
endif()

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

add_executable(dspecm1d_gia3d_export_isotropic_prem
    gia3d/export_isotropic_prem.f90
)
target_link_libraries(dspecm1d_gia3d_export_isotropic_prem PRIVATE
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

set(isotropic_prem_model_dir
    "${CMAKE_CURRENT_BINARY_DIR}/isotropic-prem-models"
)
set(isotropic_prem_models)
foreach(spacing IN ITEMS
        "100km:100000"
        "50km:50000"
        "25km:25000"
        "12.5km:12500")
    string(REPLACE ":" ";" spacing_fields "${spacing}")
    list(GET spacing_fields 0 spacing_label)
    list(GET spacing_fields 1 spacing_metres)
    set(dspec_model
        "${isotropic_prem_model_dir}/isotropic_prem_${spacing_label}.dspec"
    )
    set(gia3d_model
        "${isotropic_prem_model_dir}/isotropic_prem_${spacing_label}.gia3d"
    )
    add_custom_command(
        OUTPUT "${dspec_model}" "${gia3d_model}"
        COMMAND
            ${CMAKE_COMMAND} -E make_directory "${isotropic_prem_model_dir}"
        COMMAND
            $<TARGET_FILE:dspecm1d_gia3d_export_isotropic_prem>
            ${spacing_metres} "${dspec_model}" "${gia3d_model}"
        DEPENDS dspecm1d_gia3d_export_isotropic_prem
        VERBATIM
    )
    list(APPEND isotropic_prem_models "${dspec_model}" "${gia3d_model}")
endforeach()
add_custom_target(dspecm1d_isotropic_prem_models ALL
    DEPENDS ${isotropic_prem_models}
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

add_test(
    NAME LoveNumbersGia3DReference.ReportsIsotropicPremComparison
    COMMAND
        ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/compare_isotropic_prem.py"
        "$<TARGET_FILE:dspecm1d_controlled_love_numbers>"
        "$<TARGET_FILE:dspecm1d_gia3d_controlled_reference>"
        "$<TARGET_FILE:dspecm1d_gia3d_export_isotropic_prem>"
        "${isotropic_prem_model_dir}"
)
set_tests_properties(
    LoveNumbersGia3DReference.ReportsIsotropicPremComparison
    PROPERTIES
        LABELS "paper-validation"
        TIMEOUT 900
)

set(elln_validator
    ${Python3_EXECUTABLE}
    "${CMAKE_CURRENT_SOURCE_DIR}/elln/validate_elln.py"
    --reference
    "${CMAKE_CURRENT_SOURCE_DIR}/elln/official_ti_prem_table_s7.txt"
    --driver "${CMAKE_CURRENT_SOURCE_DIR}/elln/run_ti_prem_example.m"
    --patch "${CMAKE_CURRENT_SOURCE_DIR}/elln/patches/non_gui.patch"
    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/elln"
)
set(elln_source_arguments)
if(DSPECM1D_ELLN_ARCHIVE)
    list(APPEND elln_source_arguments
        --archive "${DSPECM1D_ELLN_ARCHIVE}"
    )
elseif(DSPECM1D_ELLN_SOURCE_DIR)
    list(APPEND elln_source_arguments
        --source-dir "${DSPECM1D_ELLN_SOURCE_DIR}"
    )
endif()
list(APPEND elln_validator ${elln_source_arguments})

find_program(DSPECM1D_MATLAB_EXECUTABLE NAMES matlab)
find_program(DSPECM1D_OCTAVE_EXECUTABLE NAMES octave octave-cli)
if(DSPECM1D_MATLAB_EXECUTABLE)
    list(APPEND elln_validator
        --matlab "${DSPECM1D_MATLAB_EXECUTABLE}"
    )
elseif(DSPECM1D_OCTAVE_EXECUTABLE)
    list(APPEND elln_validator
        --octave "${DSPECM1D_OCTAVE_EXECUTABLE}"
    )
endif()

add_custom_target(dspecm1d_validate_elln
    COMMAND ${elln_validator}
    COMMENT "Auditing the official ELLN reference"
    VERBATIM
)

add_test(
    NAME LoveNumbersELLNReference.ValidatesConversion
    COMMAND
        ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/elln/validate_elln.py"
        --reference
        "${CMAKE_CURRENT_SOURCE_DIR}/elln/official_ti_prem_table_s7.txt"
        --conversion-only
)
set_tests_properties(
    LoveNumbersELLNReference.ValidatesConversion
    PROPERTIES LABELS "paper-validation"
)

if(DSPECM1D_ELLN_ARCHIVE OR DSPECM1D_ELLN_SOURCE_DIR)
    set(elln_gravity_copy
        "${CMAKE_CURRENT_BINARY_DIR}/elln-gravity-copy"
    )
    file(MAKE_DIRECTORY
        "${elln_gravity_copy}/src"
        "${elln_gravity_copy}/validation"
    )
    foreach(source IN ITEMS
            LoveDofMap.h
            LoveNumbers.cpp
            RadialState.h
            StaticForcing.h
            StaticOperator.h
            StaticSolve.h)
        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/../src/${source}"
            "${elln_gravity_copy}/src/${source}"
            COPYONLY
        )
    endforeach()
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/controlled_love_numbers.cpp"
        "${elln_gravity_copy}/validation/controlled_love_numbers.cpp"
        COPYONLY
    )
    execute_process(
        COMMAND
            ${PATCH_EXECUTABLE} --batch --forward
            --no-backup-if-mismatch -p1
            --input=${CMAKE_CURRENT_SOURCE_DIR}/elln/patches/elln_gravitational_constant.patch
        WORKING_DIRECTORY "${elln_gravity_copy}"
        RESULT_VARIABLE elln_gravity_patch_result
        OUTPUT_VARIABLE elln_gravity_patch_output
        ERROR_VARIABLE elln_gravity_patch_error
    )
    if(NOT elln_gravity_patch_result EQUAL 0)
        message(FATAL_ERROR
            "Could not patch the validation Love-number copy:\n"
            "${elln_gravity_patch_output}${elln_gravity_patch_error}"
        )
    endif()

    add_executable(dspecm1d_controlled_love_numbers_elln_g
        "${elln_gravity_copy}/validation/controlled_love_numbers.cpp"
    )
    target_include_directories(
        dspecm1d_controlled_love_numbers_elln_g PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/../include"
    )
    target_link_libraries(
        dspecm1d_controlled_love_numbers_elln_g PRIVATE DSpecM1D
    )

    set(elln_ti_model_dir
        "${CMAKE_CURRENT_BINARY_DIR}/elln-ti-prem-models"
    )
    set(elln_ti_models
        "${elln_ti_model_dir}/elln_ti_prem_100km.dspec"
        "${elln_ti_model_dir}/elln_ti_prem_50km.dspec"
        "${elln_ti_model_dir}/elln_ti_prem_25km.dspec"
        "${elln_ti_model_dir}/elln_ti_prem_12.5km.dspec"
        "${elln_ti_model_dir}/elln_ti_prem_6.25km.dspec"
    )
    set(elln_ti_manifest
        "${elln_ti_model_dir}/elln_ti_prem_manifest.json"
    )
    add_custom_command(
        OUTPUT ${elln_ti_models} "${elln_ti_manifest}"
        COMMAND
            ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/export_ti_prem.py"
            ${elln_source_arguments}
            --output-dir "${elln_ti_model_dir}"
            --expected-checksums
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/ti_prem_models.sha256"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/export_ti_prem.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/validate_elln.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/ti_prem_models.sha256"
        VERBATIM
    )
    add_custom_target(dspecm1d_elln_ti_prem_models ALL
        DEPENDS ${elln_ti_models} "${elln_ti_manifest}"
    )

    add_test(
        NAME LoveNumbersELLNReference.AuditsOfficialExample
        COMMAND ${elln_validator}
    )
    set_tests_properties(
        LoveNumbersELLNReference.AuditsOfficialExample
        PROPERTIES
            LABELS "paper-validation"
            TIMEOUT 300
    )
    add_test(
        NAME LoveNumbersELLNReference.AuditsDegreeOneConvention
        COMMAND
            ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/validate_elln.py"
            --reference
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/official_ti_prem_table_s7.txt"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/elln-degree-one"
            --degree-one-only
            ${elln_source_arguments}
    )
    set_tests_properties(
        LoveNumbersELLNReference.AuditsDegreeOneConvention
        PROPERTIES LABELS "paper-validation"
    )
    add_test(
        NAME LoveNumbersELLNReference.ReportsTiPremComparison
        COMMAND
            ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/compare_ti_prem.py"
            "$<TARGET_FILE:dspecm1d_controlled_love_numbers>"
            "${elln_ti_model_dir}"
            "${elln_ti_manifest}"
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/ti_prem_models.sha256"
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/official_ti_prem_table_s7.txt"
    )
    set_tests_properties(
        LoveNumbersELLNReference.ReportsTiPremComparison
        PROPERTIES
            LABELS "paper-validation"
            TIMEOUT 900
    )
    add_test(
        NAME LoveNumbersELLNReference.ReportsGravitationalConstantSensitivity
        COMMAND
            ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/compare_gravitational_constants.py"
            "$<TARGET_FILE:dspecm1d_controlled_love_numbers>"
            "$<TARGET_FILE:dspecm1d_controlled_love_numbers_elln_g>"
            "${elln_ti_model_dir}/elln_ti_prem_6.25km.dspec"
            "${elln_ti_manifest}"
            "${CMAKE_CURRENT_SOURCE_DIR}/elln/official_ti_prem_table_s7.txt"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/elln-gravity-audit"
            ${elln_source_arguments}
    )
    set_tests_properties(
        LoveNumbersELLNReference.ReportsGravitationalConstantSensitivity
        PROPERTIES
            LABELS "paper-validation"
            TIMEOUT 900
    )
else()
    message(STATUS
        "ELLN source not configured: set DSPECM1D_ELLN_ARCHIVE or "
        "DSPECM1D_ELLN_SOURCE_DIR to enable official-example reproduction"
    )
endif()
