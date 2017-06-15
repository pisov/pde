#!/bin/bash

rm -f *.png

FILES=`ls -1 wave-??????.dat`
EFILE="energy.dat"


for file in $FILES
do
   printf 'Processing %s ...\r' $file
   ID=`echo $file |  sed 's/wave-//' | sed 's/.dat//'`
   fnamepng=wave-$ID.png
   echo "set size square" > frm.gnu
   echo "set yrange [-1.1:1.1]" >> frm.gnu
   echo "set xrange [0:1]" >> frm.gnu 
   printf 'set xlabel "x"\n' >> frm.gnu
   printf 'set  ylabel "u(x, t)"\n' >> frm.gnu
   time=`fgrep -e"#" $file | head -1 | awk '{print $2}'`
   printf 'set title "Time %.1f"\n' $time >> frm.gnu 
   echo "set term png" >> frm.gnu
   printf 'set output "%s"\n' $fnamepng >> frm.gnu
   printf 'plot "%s" u 1:2 w l  lw 2 notitle\n' $file >> frm.gnu
   gnuplot frm.gnu
done

printf '\nCreate animated GIF wave.gif ...'
convert *.png wave.gif
printf 'Done !!!\n'
