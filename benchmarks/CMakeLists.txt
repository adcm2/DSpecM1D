add_executable(bench_lf bench_lf.cpp)
target_link_libraries(bench_lf DSpecM1D)

add_executable(bench_mf bench_mf.cpp)
target_link_libraries(bench_mf DSpecM1D)

add_executable(bench_bolivia bench_bolivia.cpp)
target_link_libraries(bench_bolivia DSpecM1D)

add_executable(bench_step bench_step.cpp)
target_link_libraries(bench_step DSpecM1D)

add_executable(tidal_forcing tidal_forcing.cpp)
target_link_libraries(tidal_forcing DSpecM1D)

add_executable(example_record example_record.cpp)
target_link_libraries(example_record DSpecM1D)