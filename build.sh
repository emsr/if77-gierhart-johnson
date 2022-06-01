#! /bin/bash

gfortran -g -fno-align-commons -o atoa src/*.for 2> warn.txt

./atoa < atoa.in > atoa.out
