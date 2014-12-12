#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $0)

rm /tmp/qprof.pdf
rbox n D10 30 t | PROFILEFREQUENCY=10000 LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=/tmp/quickhull.profile $PROJECT_ROOT/bin/quickhull > /dev/null
google-pprof --pdf $PROJECT_ROOT/bin/quickhull /tmp/quickhull.profile > /tmp/qprof.pdf
