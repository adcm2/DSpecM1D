add_executable(testexample testexample.cpp)
target_link_libraries(testexample new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

add_executable(convergence_number convergence_number.cpp)
target_link_libraries(convergence_number new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

add_executable(augment_test augment_test.cpp)
target_link_libraries(augment_test new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

add_executable(extractmineos extractmineos.cpp)
target_link_libraries(extractmineos new_coupling PlanetaryModel FFTWpp TomographyModels)

add_executable(splinetest splinetest.cpp)
target_link_libraries(splinetest new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

add_executable(onedexample onedexample.cpp)
target_link_libraries(onedexample new_coupling PlanetaryModel FFTWpp TomographyModels gplspec Spectra)

add_executable(augment_basis augment_basis.cpp)
target_link_libraries(augment_basis new_coupling PlanetaryModel FFTWpp TomographyModels gplspec Spectra)

add_executable(toroidal toroidal.cpp)
target_link_libraries(toroidal new_coupling)
# # target_precompile_headers(toroidal REUSE_FROM new_coupling)
# target_precompile_headers(toroidal PRIVATE "pch.h")

add_executable(toroidal_bench toroidal_bench.cpp)
target_link_libraries(toroidal_bench new_coupling)

add_executable(toroidal_sphere toroidal_sphere.cpp)
target_link_libraries(toroidal_sphere new_coupling)

add_executable(toroidal_mineos toroidal_mineos.cpp)
target_link_libraries(toroidal_mineos new_coupling)
# target_precompile_headers(toroidal_mineos PRIVATE "pch.h")

add_executable(toroidal_gal toroidal_gal.cpp)
target_link_libraries(toroidal_gal new_coupling)
# target_precompile_headers(toroidal_gal PRIVATE "pch.h")

add_executable(toroidal_spectra toroidal_spectra.cpp)
target_link_libraries(toroidal_spectra new_coupling GSHTrans)
# target_precompile_headers(toroidal_spectra PRIVATE "pch.h")
# Add this line after add_executable() or add_library()

# target_precompile_headers(toroidal REUSE_FROM new_coupling)
# target_precompile_headers(toroidal PRIVATE "pch.h")

add_executable(toroidal_sem toroidal_sem.cpp)
target_link_libraries(toroidal_sem new_coupling GSHTrans)
# target_precompile_headers(toroidal_sem PRIVATE pch2.h)
# target_compile_options(toroidal_sem PRIVATE -ftime-report)

add_executable(test_spectra test_spectra.cpp)
target_link_libraries(test_spectra new_coupling PlanetaryModel FFTWpp TomographyModels gplspec Spectra)

add_executable(extracttoroidal extracttoroidal.cpp)
target_link_libraries(extracttoroidal new_coupling PlanetaryModel FFTWpp TomographyModels)

add_executable(premrewrite premrewrite.cpp)
target_link_libraries(premrewrite new_coupling PlanetaryModel FFTWpp TomographyModels)

add_executable(testwigner testwigner.cpp)
target_link_libraries(testwigner new_coupling PlanetaryModel FFTWpp TomographyModels)

add_executable(lsem_toroidal lsem_toroidal.cpp)
target_link_libraries(lsem_toroidal new_coupling GSHTrans)

add_executable(toroidal_l_augtest toroidal_l_augtest.cpp)
target_link_libraries(toroidal_l_augtest new_coupling GSHTrans)

add_executable(toroidal_comparison toroidal_comparison.cpp)
target_link_libraries(toroidal_comparison new_coupling GSHTrans)