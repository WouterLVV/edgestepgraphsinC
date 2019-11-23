# cmake_minimum_required(VERSION <specify CMake version here>)
project(edgestepgraph C)

set(CMAKE_C_STANDARD 11)

add_executable(edgestepgraph main.c)


add_executable(edgestepgraph2 alt.c)
target_link_libraries(edgestepgraph2 m)


add_executable(edgestepgraph3 alt2.c)
target_link_libraries(edgestepgraph3 m)

add_executable(edgestepgraph4 alt3.c)
target_link_libraries(edgestepgraph4 m)