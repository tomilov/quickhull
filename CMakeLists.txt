project("quickhull")

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories("include/")

add_executable("${PROJECT_NAME}" "src/quickhull.cpp"  "include/quickhull.hpp")
add_executable("qh"              "src/simple_use.cpp" "include/quickhull.hpp")
