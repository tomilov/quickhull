project("quickhull")

cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

add_compile_options(-W)
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    add_compile_options(-Weverything)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wsuggest-override)
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER_ID} currently is not supported")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG=1 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -fno-inline -DDEBUG=1 -march=x86-64 -mtune=generic")
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -D_GLIBCXX_DEBUG -D_LIBCPP_DEBUG2=1")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
endif()
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -DNDEBUG=1 -march=native")
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -gline-tables-only")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Og -ggdb")
endif()
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os -DNDEBUG=1 -march=native")

include_directories("include/")

add_executable("${PROJECT_NAME}" "src/quickhull.cpp"  "include/quickhull.hpp")
#add_executable("qh"              "src/simple_use.cpp" "include/quickhull.hpp")

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
    )
