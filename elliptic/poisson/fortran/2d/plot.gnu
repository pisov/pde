set xlabel 'x'
set ylabel 'y'
set dgrid3d 100, 100
set hidden3d
splot 'plot.dat' u 1:2:3 notitle w l
pause -1
