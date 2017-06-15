#!/bin/bash

rm -f *.png

START=0
END=50000
STEP=100

for i in `seq $START $STEP $END`
do
   fname=`printf 'wf-%.5d.dat' $i`
   fnamepng=`printf 'wf-%.5d.png' $i`
   echo "set size square" > frm.gnu
   echo "set yrange [0:0.1]" >> frm.gnu
   echo "set xrange [-1:1]" >> frm.gnu 
   echo "set term png" >> frm.gnu
   printf 'set output "%s"\n' $fnamepng >> frm.gnu
   printf 'plot "%s" u 1:($2**2+$3**2) w l notitle\n' $fname >> frm.gnu
   gnuplot frm.gnu
done

convert -delay 50 *.png wf.gif

