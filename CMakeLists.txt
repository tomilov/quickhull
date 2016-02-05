project(quickhull)

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)

if(NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(FATAL_ERROR "only clang supported currently")
endif()

include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" COMPILER_SUPPORTS_LIBCXX)
if(NOT COMPILER_SUPPORTS_LIBCXX)
    message(FATAL_ERROR "libc++ not installed or libc++ not supported by clang++")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")

CHECK_CXX_COMPILER_FLAG("-std=gnu++1z" COMPILER_SUPPORTS_CXX1Z)
if(NOT COMPILER_SUPPORTS_CXX1Z)
    message(FATAL_ERROR "Compiler does not support C++1z standard")
endif()
add_compile_options("-std=gnu++1z")

#set(CXX_EXTENSIONS ON)
#set_property(TARGET ${TARGETS} PROPERTY CXX_STANDARD 17)
#set_property(TARGET ${TARGETS} PROPERTY CXX_STANDARD_REQUIRED ON)
add_compile_options(-fconstexpr-steps=50000000 -ftemplate-depth=200)
add_compile_options(-ftemplate-backtrace-limit=0 -fconstexpr-backtrace-limit=0 -fdiagnostics-show-template-tree)
add_compile_options(-W -Weverything -Wmissing-include-dirs -Wconversion -Wmismatched-tags -Wuninitialized -Wfatal-errors -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded)
add_compile_options(-Wno-gnu-anonymous-struct -Wno-nested-anon-types)
#add_compile_options(-ftime-report)

set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -fno-inline -DDEBUG=1 -D_DEBUG=1 -D_GLIBCXX_DEBUG=1 -march=x86-64 -mtune=generic")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG=1 -march=native")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -gline-tables-only -DNDEBUG=1 -march=native")
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG=1 -march=native")
else()
    message(STATUS "Wrong build type selected, default to Debug.")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "") # suppress -rdynamic flag

include_directories("include/")

add_executable("${PROJECT_NAME}" "src/quickhull.cpp"  "include/quickhull.hpp")
add_executable("qh"              "src/simple_use.cpp" "include/quickhull.hpp")

set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
install(TARGETS "${PROJECT_NAME}" "qh" RUNTIME DESTINATION "bin/")
