# add_executable(spheroidal_1 spheroidal_1.cpp)
# target_link_libraries(spheroidal_1 new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

# add_executable(spheroidal_lmax spheroidal_lmax.cpp)
# target_link_libraries(spheroidal_lmax new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

# add_executable(tensorcalc tensorcalc.cpp)

# add_executable(full_lmax full_lmax.cpp)
# target_link_libraries(full_lmax new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

add_executable(debug_sph debug_sph.cpp)
target_link_libraries(debug_sph new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

# add_executable(debug_tor debug_tor.cpp)
# target_link_libraries(debug_tor new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

# add_executable(debug_rad debug_rad.cpp )
# target_link_libraries(debug_rad new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

# add_executable(debug_full debug_full.cpp)
# target_link_libraries(debug_full new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)   

# add_executable(debug_processing debug_processing.cpp)
# target_link_libraries(debug_processing new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

add_executable(debug_atten debug_atten.cpp)
target_link_libraries(debug_atten new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

# add_executable(bench_highfrequency bench_highfrequency.cpp)
# target_link_libraries(bench_highfrequency new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

# add_executable(bench_lf bench_lf.cpp)
# target_link_libraries(bench_lf new_coupling PlanetaryModel FFTWpp TomographyModels gplspec OpenMP::OpenMP_CXX)

add_executable(paper_examples paper_examples.cpp)
target_link_libraries(paper_examples new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)

add_executable(example_record example_record.cpp)
target_link_libraries(example_record new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)
# add_executable(test_start test_start.cpp)
# target_link_libraries(test_start new_coupling PlanetaryModel FFTWpp TomographyModels gplspec)
add_executable(check_script check_script.cpp)
target_link_libraries(check_script new_coupling PlanetaryModel)
