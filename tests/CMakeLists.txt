# tests/CMakeLists.txt

# 1. Fetch GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0 # Latest stable version
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# 2. Create your test executable
add_executable(dspecm1d_tests 
    test_math.cpp
    # Add more test files here as your project grows: test_mesh.cpp, test_fft.cpp, etc.
)

# 3. Link GoogleTest AND your main library (DSpecM1D)
target_link_libraries(dspecm1d_tests PRIVATE 
    DSpecM1D 
    gtest_main # gtest_main automatically provides the main() function!
)

# 4. Tell CMake to discover and register your GoogleTests automatically
include(GoogleTest)
gtest_discover_tests(dspecm1d_tests)