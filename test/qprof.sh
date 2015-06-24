#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $0)/..

rm /tmp/qprof.pdf
rbox n D10 38 s t | PROFILEFREQUENCY=10000 LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=/tmp/quickhull.profile /tmp/quickhull
google-pprof --pdf /tmp/quickhull /tmp/quickhull.profile > /tmp/qprof.pdf
gnome-open /tmp/qprof.pdf
