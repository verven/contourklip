cmake_minimum_required(VERSION 3.20)
project(contourklip)

set(CMAKE_CXX_STANDARD 20)
set(BINARY ${CMAKE_PROJECT_NAME})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer")

include_directories(include)
include_directories(single_include)
include_directories(test/thirdparty)

add_executable(${BINARY}_unittst test/unit.cpp
        test/test_utilities.hpp
#        ${include_files}
        include/geometry_base.hpp
        include/bezier_utils.hpp
        include/svg_io.hpp
        include/sweep_point.hpp
        include/polyclip.hpp
        )
add_executable(${BINARY}_acceptancetst test/acceptance.cpp)
add_executable(${BINARY}_test_svg test/svg_testcases.cpp)
add_executable(${BINARY}_example examples/example1.cpp)
add_executable(${BINARY}_example2 examples/example2.cpp)