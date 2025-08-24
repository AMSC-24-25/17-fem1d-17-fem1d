cmake_minimum_required(VERSION 3.26)
project(17_fem1d_17_fem1d)

set(CMAKE_CXX_STANDARD 20)

# ---- GoogleTest-----
add_subdirectory(googletest)
enable_testing()


# --------include ---------

include_directories(include)
include_directories(include/external)  # For toml.hpp and exprtk.hpp

set(CMAKE_PREFIX_PATH "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/share/eigen3/cmake/")
find_package(Eigen3 3.3 REQUIRED)
include_directories(${mkEigenInc})

# OpenMP support (optional)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - enabling parallel support")
else()
    message(WARNING "OpenMP not found - building without parallel support")
endif()

# ----- Sources -----

file(GLOB ALL_SRC_FILES src/*.cpp)
list(FILTER ALL_SRC_FILES EXCLUDE REGEX ".*config.cpp$")
list(FILTER ALL_SRC_FILES EXCLUDE REGEX ".*TomlMain.cpp$")
list(FILTER ALL_SRC_FILES EXCLUDE REGEX ".*main.cpp$")
list(FILTER ALL_SRC_FILES EXCLUDE REGEX ".*main_td.cpp$")

# Create external library list
set(FEM_LIBRARIES Eigen3::Eigen)

if(OpenMP_CXX_FOUND)
    message("OpenMP was found. Compiling sequential library as well as parallel library.")
    add_library(fem1d_lib_sequential STATIC ${ALL_SRC_FILES})
    target_link_libraries(fem1d_lib_sequential ${FEM_LIBRARIES}) 
    add_library(fem1d_lib STATIC ${ALL_SRC_FILES})
    target_link_libraries(fem1d_lib ${FEM_LIBRARIES} OpenMP::OpenMP_CXX)
else()
    message("OpenMP was not found. Compiling library as sequential.")
    add_library(fem1d_lib STATIC ${ALL_SRC_FILES})
    target_link_libraries(fem1d_lib ${FEM_LIBRARIES})
endif()


# Eseguibile principale
add_executable(17_fem1d_17_fem1d src/main.cpp)
target_link_libraries(17_fem1d_17_fem1d fem1d_lib)

# Nuovo eseguibile con configurazione TOML
add_executable(TomlMain src/TomlMain.cpp src/config.cpp)
target_link_libraries(TomlMain fem1d_lib)

if(OpenMP_CXX_FOUND)
    message("Compiling sequential TomlMain for speedup analysis.")
    add_executable(sequentialTomlMain src/TomlMain.cpp src/config.cpp)
    target_link_libraries(sequentialTomlMain fem1d_lib_sequential)
endif()

# Eseguibile per main_td.cpp
add_executable(17_fem1d_17_fem1d_td src/main_td.cpp)
target_link_libraries(17_fem1d_17_fem1d_td fem1d_lib)

# ---- Output preparation ----
add_custom_target(create_output_dir ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory output
) 
add_custom_target(empty_output_dir ALL 
        COMMAND ${CMAKE_COMMAND} -E remove
        "output/*.vtu"
        "output/*.csv"
        COMMENT "Removing all .vtu and .csv files from the output directory" 
) 
add_dependencies(17_fem1d_17_fem1d create_output_dir)
add_dependencies(17_fem1d_17_fem1d_td create_output_dir)
add_dependencies(TomlMain create_output_dir)
if(OpenMP_CXX_FOUND)
    add_dependencies(sequentialTomlMain create_output_dir)
endif()

# -------testing----

# Test - Auto-include all test files
file(GLOB TEST_SRC_FILES test/*.cpp)
add_executable(17_fem1d_17_fem1d_test ${TEST_SRC_FILES})
# Collegamento alla libreria invece di includere i file .cpp direttamente
target_link_libraries(17_fem1d_17_fem1d_test fem1d_lib gtest_main)

include(GoogleTest)
gtest_discover_tests(17_fem1d_17_fem1d_test)
