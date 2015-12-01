cmake_minimum_required(VERSION 3.2)
project(mytool)


add_subdirectory(librarymsms)
add_subdirectory(groundtruth)
add_subdirectory(calcrecallprecision)
add_subdirectory(readpepxml)
add_subdirectory(src)
add_subdirectory(swathnoise)
add_subdirectory(RemoveNoise)

