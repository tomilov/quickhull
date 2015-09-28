project( quickhull )
cmake_minimum_required( VERSION 2.8 )

if(     CMAKE_BUILD_TYPE STREQUAL "Debug" )
elseif( CMAKE_BUILD_TYPE STREQUAL "Release" )
elseif( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
elseif( CMAKE_BUILD_TYPE STREQUAL "MinSizeRel" )
else()
    message( STATUS "No build type selected, default to MinSizeRel" )
    set( CMAKE_BUILD_TYPE "MinSizeRel" )
endif()

set( CMAKE_VERBOSE_MAKEFILE ON )
set( CMAKE_COLOR_MAKEFILE   ON )

set( CMAKE_CXX_COMPILER             "clang++"                                                                 )
set( CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -W -Weverything -pedantic -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -ftemplate-backtrace-limit=0 -fdiagnostics-color=always" )
set( CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -std=gnu++1z -stdlib=libc++"                          )
set( CMAKE_CXX_FLAGS_RELEASE        "-Ofast -DNDEBUG=1 -march=native"                                         )
set( CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG=1 -march=native"                                            )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -gline-tables-only -DNDEBUG=1 -march=native"                         )
set( CMAKE_CXX_FLAGS_DEBUG          "-g -DDEBUG=1 -D_DEBUG=1 -D_GLIBCXX_DEBUG=1 -march=x86-64 -mtune=generic" )

set( CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "" ) # suppress -rdynamic flag

#aux_source_directory( src/ SOURCES )
include_directories( include/ )

add_executable( ${PROJECT_NAME} "src/quickhull.cpp" "include/quickhull.hpp" )
