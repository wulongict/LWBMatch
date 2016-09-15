cmake_minimum_required(VERSION 3.2)
project(libmsms)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fpermissive -fopenmp ")
set (CMAKE_CXX_LINK_FLAGS "-lz")


set(SOURCE_FILES
        PeakList.cpp
        SpectraST_base64.cpp
        SpectraST_cramp.cpp
        SpectraST_ramp.cpp
        Util.cpp
        SWATHmzXMLReader.cpp
        mzXMLReader.cpp
        Util.h
)




add_library(msms SHARED ${SOURCE_FILES})
target_link_libraries(msms -lz -pthread)

install(TARGETS msms DESTINATION ${PROJECT_SOURCE_DIR}/../bin)
