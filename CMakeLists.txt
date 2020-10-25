cmake_minimum_required(VERSION 3.17)
project(pazi_lab_1 C)

include_directories(/usr/local/Cellar/libgcrypt/1.8.6/include)
link_directories(/usr/local/Cellar/libgcrypt/1.8.6/lib)
link_libraries(gcrypt)

set(CMAKE_C_STANDARD 11)

add_executable(pazi_lab_1 main.c jacobi_curve/jacobi_curve.h)