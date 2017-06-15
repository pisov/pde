1. Compile

gfortran -O3 poisson_2d.f90 -o poisson_2d.x

2. Execute the code

./poisson_2d.x > plot.dat

n = 100
m = 100

3. Plot the result data

gnuplot plot.gnu 
