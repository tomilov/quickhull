#!/usr/bin/env bash

rbox n D10 30 s > /tmp/quickhull.txt
LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=/tmp/quickhull.profile bin/quickhull /tmp/quickhull.txt
google-pprof --pdf bin/quickhull /tmp/quickhull.gprof > qprof.pdf
