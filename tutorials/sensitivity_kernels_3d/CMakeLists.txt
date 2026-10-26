cmake_minimum_required(VERSION 3.24)
project(DSpecM1D3DKernelConsumer LANGUAGES CXX)
find_package(DSpecM1D CONFIG REQUIRED)
add_executable(dspec_3d_kernel_consumer consumer.cpp)
target_link_libraries(dspec_3d_kernel_consumer PRIVATE DSpecM1D::DSpecM1D)

add_executable(dspec_3d_kernel_timing kernel_timing.cpp)
target_link_libraries(dspec_3d_kernel_timing PRIVATE DSpecM1D::DSpecM1D)
