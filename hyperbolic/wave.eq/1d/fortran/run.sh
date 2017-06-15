#!/bin/bash

rm -f *.x && gfortran -O3 wave.f90 -o wave.x && rm -f *.dat && ./wave.x < input
