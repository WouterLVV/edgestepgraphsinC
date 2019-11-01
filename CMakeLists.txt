# cmake_minimum_required(VERSION <specify CMake version here>)
project(edgestepgraph C)

set(CMAKE_C_STANDARD 11)

add_executable(edgestepgraph main.c)
add_executable(edgestepgraph2 alt.c)