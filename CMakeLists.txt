cmake_minimum_required(VERSION 3.2)
project(expfloat)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(include)

set(HEADER_FILES include/expansion_math.h include/utils.h include/test_apps.h)
set(SOURCE_FILES src/main.cpp)

add_executable(expfloat ${SOURCE_FILES} ${HEADER_FILES})