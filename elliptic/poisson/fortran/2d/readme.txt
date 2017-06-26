1. Compile

gfortran -O3 poisson.f90 -o poisson.x

2. Execute the code

./poisson.x > plot.dat

n = 100
m = 100

3. Plot the result data

gnuplot plot.gnu 
