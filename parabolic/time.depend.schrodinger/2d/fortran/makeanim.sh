#!/bin/bash

rm -f *.png


FILES=`ls -1 wf-??????.ppm`

convert potEnergy.ppm -transparent blue potEnergy.png
convert  potEnergy.png -fill black -opaque red potEnergyBlack.png

for file in $FILES
do
   ID=`echo $file |  sed 's/wf-//' | sed 's/.ppm//'`
   fnamepng=wf-$ID.png
   composite -gravity center potEnergyBlack.png $file $fnamepng   
done

rm -f $$

convert -delay 10 wf-*.png wf.gif
ffmpeg -i wf.gif -strict -2  -pix_fmt yuv420p -y wf.mp4
