#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $0)

rm /tmp/qprof.pdf
rbox n D10 30 s > /tmp/quickhull.txt
LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=/tmp/quickhull.profile $PROJECT_ROOT/bin/quickhull /tmp/quickhull.txt
google-pprof --pdf $PROJECT_ROOT/bin/quickhull /tmp/quickhull.profile > /tmp/qprof.pdf
