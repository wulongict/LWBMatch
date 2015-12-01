cmake_minimum_required(VERSION 3.2)
project(groundtruth)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    ExtractGroundTruth.cpp
    ../librarymsms/FileReader.h
    )

add_executable(groundtruth ${SOURCE_FILES})