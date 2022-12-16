cmake_minimum_required(VERSION 3.22)
project(MyVector)

set(CMAKE_CXX_STANDARD 20)

find_package(MPI REQUIRED)

add_executable(MyVector main.cpp Matrix.cpp Matrix.h )
target_link_libraries(MyVector PUBLIC MPI::MPI_CXX)
add_executable(Generator generator.cpp)
