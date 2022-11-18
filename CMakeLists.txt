cmake_minimum_required(VERSION 3.10)
project(pcog)

set(COMPILE_FLAGS "-Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Wcast-align -Wold-style-cast
 -Wunused -Woverloaded-virtual -Wconversion -Wpedantic -Wmisleading-indentation
-Wimplicit-fallthrough -Wformat=2 -Wdouble-promotion -Wuseless-cast -Wnull-dereference -Wlogical-op -Wduplicated-cond -Wduplicated-branches")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILE_FLAGS}")
