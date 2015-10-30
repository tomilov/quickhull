#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $0)/..

BUILD_DIR=$PROJECT_ROOT/bin/
rm /tmp/qprof.pdf
rbox n D10 35 s t | PROFILEFREQUENCY=10000 LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=$BUILD_DIR/quickhull.profile $BUILD_DIR/quickhull
google-pprof --pdf $BUILD_DIR/quickhull $BUILD_DIR/quickhull.profile > $BUILD_DIR/qprof.pdf
gnome-open $BUILD_DIR/qprof.pdf
