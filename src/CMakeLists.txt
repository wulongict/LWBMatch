cmake_minimum_required(VERSION 3.2)
project(lwbmatch)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fpermissive -Wunused-result -fopenmp ")
set (CMAKE_CXX_LINK_FLAGS "-lz")

set(SOURCE_FILES
    alignmentevaluation.cpp
    alignmentevaluation.h
    dtw.cpp
    dtw.h
    featuremap.cpp
    featuremap.h
    lowess.cpp
    lowess.h
    lwbmatch.cpp
    lwbmatch.h
    main.cpp
    SparseMatrix.h)

#set( PROJECT_LINK_LIBS libmsms.so )
link_directories( ../librarymsms/ ) # lib1.so and lib2.so are there.
link_directories(/usr/local/lib/)

add_executable(lwbmatch ${SOURCE_FILES} Dewarper.cpp Dewarper.h parameters.cpp parameters.h)
target_link_libraries(lwbmatch -lmsms)
add_dependencies(lwbmatch msms)



add_executable(mzXMLReader ReadmzXML.cpp)
target_link_libraries(mzXMLReader -lz -pthread -lmsms -lboost_program_options)
#target_link_libraries(lwbmatch -lmsms)

add_dependencies(mzXMLReader msms)


add_executable(mzxmltomatrix mzXMLtomatrix.cpp)
target_link_libraries(mzxmltomatrix -lz -pthread -lmsms -lboost_program_options)

add_dependencies(mzxmltomatrix msms)
