cmake_minimum_required(VERSION 3.2)
project(calcrecallprecision)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    CalculateRecallPrecision.cpp)

add_executable(recallprecision ${SOURCE_FILES})
install(TARGETS recallprecision DESTINATION ${PROJECT_SOURCE_DIR}/../bin)
