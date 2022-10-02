cmake_minimum_required(VERSION 3.21)
project (MPI)

set(CMAKE_CXX_STANDARD 14)

# Add MPI Package to Project
find_package(MPI REQUIRED)

add_executable(MPI main.cpp)
# Add libraries for code completion and compiling
target_link_libraries(MPI PUBLIC MPI::MPI_CXX)
