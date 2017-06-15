#!/bin/bash

rm -f *.png

START=0
END=50000
STEP=100

FILES=`ls -1 wf-??????.dat`
EFILE="energy.dat"


for file in $FILES
do
   ID=`echo $file |  sed 's/wf-//' | sed 's/.dat//'`
   fnamepng=wf-$ID.png
   ENERG=`tail -1 $file | awk '{printf("%f %f\n",$6,$7)}'`
   echo "set size square" > frm.gnu
   echo "set yrange [-1.1:1.1]" >> frm.gnu
   echo "set y2range [:20]" >> frm.gnu
   echo "set xrange [-25:25]" >> frm.gnu 
   echo "set ytics nomirror" >> frm.gnu
   echo "set y2tics" >> frm.gnu
   printf 'set xlabel "x"\n' >> frm.gnu
   printf 'set  ylabel "Psi(x)"\n' >> frm.gnu
   printf 'set y2label "Energy"\n' >> frm.gnu
   time=`fgrep -e"#" $file | head -1 | awk '{print $2}'`
   printf 'set title "Time %.1f"\n' $time >> frm.gnu 
   echo $time $ENERG > $$
   echo "set term png" >> frm.gnu
   printf 'set output "%s"\n' $fnamepng >> frm.gnu
   printf 'plot "%s" u 1:($4) w l  lw 2 notitle, "%s" u 1:($2) w l title "Real", "%s" u 1:($3) w l title "Imaginary", "%s" u 1:($5) w l notitle axes x1y2, "%s" u 1:2 w l lw 2 title "Kin" axes x2y2, "%s" u 1:($2+2*$3) w l lw 2 title "Pot" axes x2y2, "%s" u 1:2 w p pt 7 ps 1.5 notitle axes x2y2 , "%s" u 1:($2+2*$3) w p pt 7 ps 1.5 notitle axes x2y2\n' $file $file $file $file $EFILE $EFILE $$ $$ >> frm.gnu
   gnuplot frm.gnu
done

rm -f $$

convert -delay 10 *.png wf.gif
ffmpeg -i wf.gif -strict -2  -pix_fmt yuv420p -y wf.mp4
