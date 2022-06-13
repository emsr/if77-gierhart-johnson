#! /bin/bash

gfortran -g -fno-align-commons -o atoa src/*.for 2> warn.txt

./atoa < atoa.in > atoa.out

gfortran -g -fno-align-commons -o run_raytrace_old run_raytrace_old.for src/amsc.for src/asorp.for

gfortran -g -fno-align-commons -o run_raytrace run_raytrace.f95

g++ -g -DQUIET=0 -o run_raytrace_new run_raytrace_new.cpp
