cmake_minimum_required(VERSION 3.5)
project(solinas64)

set(CMAKE_CXX_STANDARD 11)

# Solinas64 library source files
set(SOLINAS64_LIB_SRCFILES
        solinas64.cpp
        solinas64.h)

add_library(solinas64 ${SOLINAS64_LIB_SRCFILES})

add_executable(tests tests/tests.cpp)
target_link_libraries(tests solinas64)

add_executable(benchmarks
	tests/benchmarks.cpp
	tests/gf256.h
	tests/gf256.cpp)
target_link_libraries(benchmarks solinas64)
