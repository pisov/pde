#!/bin/bash

rm -f wf* && gfortran -O3 tdse.f90 -o tdse.x -llapack && ./tdse.x < input > energy.dat && ./makeanim.sh
