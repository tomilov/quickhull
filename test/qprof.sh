#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $0)/..

rm /tmp/qprof.pdf
rbox n D10 30 s t | PROFILEFREQUENCY=10000 LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=$PROJECT_ROOT/bin/quickhull.profile ../bin/quickhull
google-pprof --pdf ../bin/quickhull $PROJECT_ROOT/bin/quickhull.profile > ../bin/qprof.pdf
gnome-open ../bin/qprof.pdf
