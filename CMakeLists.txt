cmake_minimum_required(VERSION 2.8)
project(expfloat)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fext-numeric-literals")

include_directories(include)
link_directories(/home/hari/Research/expfloat/lib)

set(HEADER_FILES include/expansion_math.h include/utils.h
  include/test_apps.h include/drecho.h)
set(SOURCE_FILES src/main.cpp src/drecho.cpp)

add_executable(expfloat ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(expfloat quadmath)

