cmake_minimum_required(VERSION 3.2)
project(LWBMatch)


add_subdirectory(librarymsms)
add_subdirectory(groundtruth)
add_subdirectory(calcrecallprecision)
add_subdirectory(src)

add_dependencies(mzxmltomatrix msms)
add_dependencies(mzXMLReader msms)
add_dependencies(lwbmatch msms)

