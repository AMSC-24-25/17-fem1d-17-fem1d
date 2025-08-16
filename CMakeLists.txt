cmake_minimum_required(VERSION 3.26)
project(17_fem1d_17_fem1d)

set(CMAKE_CXX_STANDARD 14)


#GoogleTest
add_subdirectory(googletest)
enable_testing()

include_directories(include)

set(CMAKE_PREFIX_PATH "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/share/eigen3/cmake/")
find_package(Eigen3 3.3 REQUIRED)
include_directories(${mkEigenInc})

add_executable(17_fem1d_17_fem1d
        include/boundary_cond.hpp
        include/fem1d.hpp
        include/function.hpp
        include/grid1D.hpp
        include/phi_function.hpp
        include/quadrature.hpp
        include/thomas.hpp
        src/fem1d.cpp
        src/function.cpp
        src/grid1D.cpp
        src/main.cpp
        src/phi_function.cpp
        src/quadrature.cpp
        src/thomas.cpp)

target_link_libraries(17_fem1d_17_fem1d Eigen3::Eigen)


# -------

add_executable(17_fem1d_17_fem1d_test test/main.cpp
        test/fem1d.cpp test/grid2d.cpp src/grid2D.cpp)

target_link_libraries(17_fem1d_17_fem1d_test Eigen3::Eigen gtest_main)
include(GoogleTest)

gtest_discover_tests(17_fem1d_17_fem1d_test)
