add_executable(bench_highfrequency bench_highfrequency.cpp)
target_link_libraries(bench_highfrequency new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

add_executable(bench_lf bench_lf.cpp)
target_link_libraries(bench_lf new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

add_executable(bench_mf bench_mf.cpp)
target_link_libraries(bench_mf new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

add_executable(bench_st bench_st.cpp)
target_link_libraries(bench_st new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

add_executable(bench_bolivia bench_bolivia.cpp)
target_link_libraries(bench_bolivia new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

add_executable(bench_tor bench_tor.cpp)
target_link_libraries(bench_tor new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

add_executable(bench_step bench_step.cpp)
target_link_libraries(bench_step new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

add_executable(bench_fast bench_fast.cpp)
target_link_libraries(bench_fast new_coupling)


add_executable(bench_red bench_red.cpp)
target_link_libraries(bench_red new_coupling)

add_executable(test_rot test_rot.cpp)
target_link_libraries(test_rot new_coupling)

add_executable(tidal_forcing tidal_forcing.cpp)
target_link_libraries(tidal_forcing new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)