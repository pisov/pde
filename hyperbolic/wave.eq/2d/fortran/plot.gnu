set dgrid3d 102,102
set ticslevel 0.                   
set hidden3d                       
set pm3d at s 
set palette rgbformulae 33,13,10
set xlabel "x"
set ylabel "y"
set zlabel "u(x,y)"
set xrange [0:1]
set yrange [0:1]
set zrange [0:1]
splot "wave-000000.dat" with pm3d notitle
pause -1
