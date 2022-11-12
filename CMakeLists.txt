cmake_minimum_required(VERSION 3.23)
project(slae)

set(CMAKE_CXX_STANDARD 20)

add_executable(slae main.cpp functions.cpp class.h class.cpp functions.h)
