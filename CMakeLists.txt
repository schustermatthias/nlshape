# cmake version to be used
cmake_minimum_required(VERSION 3.10)
# project name
project(nonlocal_assembly)
# Look at
# https://martinopilia.com/posts/2018/09/15/building-python-extension.html
# for Python with CMake...

# flags ----------------------------------------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)
add_compile_options("-O2")
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # Equivalent to -fPIC option in the compiler

# include --------------------------------------------------------------------------------------------------------------
include_directories(./include/ ./src/)

# target ---------------------------------------------------------------------------------------------------------------
# Adds the library Cassemble. Leads to: libCassemble.a
add_library(Cassemble SHARED
        ./src/Cassemble.cpp ./include/Cassemble.h
        ./src/MeshTypes.cpp
        ./src/mathhelpers.cpp
        ./src/model.cpp
        ./src/integration.cpp
        )

target_compile_options(Cassemble PUBLIC "-fopenmp" "-Wall")
set_target_properties(Cassemble PROPERTIES PUBLIC_HEADER "include/MeshTypes.h;include/Cassemble.h")

# external libraries ---------------------------------------------------------------------------------------------------
target_link_libraries(Cassemble PRIVATE lapack blas armadillo)
# Find OpneMP Library
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Cassemble PRIVATE OpenMP::OpenMP_CXX)
endif()
find_package(CGAL)
target_link_libraries(Cassemble PRIVATE CGAL::CGAL)

# others ---------------------------------------------------------------------------------------------------------------
add_subdirectory(./examples/CompareCodes)
# Add tests to project
option(PACKAGE_TESTS "Build the tests" ON)
#if(PACKAGE_TESTS)
#    enable_testing()
#    add_subdirectory(tests)
#endif()

# install --------------------------------------------------------------------------------------------------------------
INSTALL(TARGETS Cassemble
        LIBRARY DESTINATION $ENV{HOME}/lib
        PUBLIC_HEADER DESTINATION $ENV{HOME}/lib/include
        )
