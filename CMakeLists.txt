cmake_minimum_required(VERSION 3.15)
project(SVector_Fedchenko_Yaroslav C)

set(CMAKE_C_STANDARD 99)

add_executable(SVector_Fedchenko_Yaroslav main.c SMatrix_tests.h src/SMatrix.c src/SMatrix.h SMatrix_tests.c)