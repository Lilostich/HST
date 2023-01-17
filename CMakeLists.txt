cmake_minimum_required(VERSION 3.22)
project(MyVector)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -o myvector.o")

find_package(MPI REQUIRED)
add_library(mv main.cpp Matrix.cpp Matrix.h)
add_executable(MyVector main.cpp Matrix.cpp Matrix.h )
target_link_libraries(MyVector PUBLIC MPI::MPI_CXX)
add_executable(Generator generator.cpp)
