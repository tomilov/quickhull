#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $0)/..

rm /tmp/qprof.pdf
rbox n D10 35 s t | PROFILEFREQUENCY=10000 LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=$PROJECT_ROOT/build/quickhull.profile ../build/quickhull
google-pprof --pdf ../build/quickhull $PROJECT_ROOT/build/quickhull.profile > ../build/qprof.pdf
gnome-open ../build/qprof.pdf
