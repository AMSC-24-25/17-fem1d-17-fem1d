cmake_minimum_required(VERSION 3.26)
project(17_fem1d_17_fem1d)

set(CMAKE_CXX_STANDARD 17)


# ---- GoogleTest-----
add_subdirectory(googletest)
enable_testing()


# --------include ---------

include_directories(include)

set(CMAKE_PREFIX_PATH "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/share/eigen3/cmake/")
find_package(Eigen3 3.3 REQUIRED)
include_directories(${mkEigenInc})

# ----- Sources -----

# Crea una libreria con i file sorgente comuni
file(GLOB LIB_SRC_FILES src/*.cpp)
add_library(fem1d_lib STATIC ${LIB_SRC_FILES})
target_link_libraries(fem1d_lib Eigen3::Eigen)

# Eseguibile principale
add_executable(17_fem1d_17_fem1d src/main.cpp)
target_link_libraries(17_fem1d_17_fem1d fem1d_lib)

# -------testing----

# Test - Auto-include all test files
file(GLOB TEST_SRC_FILES test/*.cpp)
add_executable(17_fem1d_17_fem1d_test ${TEST_SRC_FILES})
# Collegamento alla libreria invece di includere i file .cpp direttamente
target_link_libraries(17_fem1d_17_fem1d_test fem1d_lib gtest_main)

include(GoogleTest)

gtest_discover_tests(17_fem1d_17_fem1d_test)
