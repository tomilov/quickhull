#!/usr/bin/env bash -vex

SRC=src/quickhull.cpp
BOOST_INCLUDE="D:/libs/boost"
INCLUDE="include/"
DEFINES="-DBOOST_VARAINT_MAX_MULTIVIZITOR_PARAMS=6"
OPT="-m64 -Ofast -march=corei7-avx -mtune=corei7-avx"

cls
time g++ $OPT -std=gnu++1y -Wall -Wextra -pedantic -isystem $BOOST_INCLUDE -I$INCLUDE $DEFINES $SRC -o bin/quickhull 2>&1
bin/quickhull
