# This example demonstrate time evolution of wave packet
#

1. Make own copy

cp tdse_2d_template.f90 tdse_2d.f90

2. Compile the source code:


gfortran -c utils.f90
gfortran -O3 tdse_2d.f90 utils.f90 -o tdse.x

or

make

3. Clear previous data files if any

rm -f wf*

4. Execute the program to evolve the wave function in time

# Exeample input values
# dx   = 0.1
# dt   = 0.0001
# kx   = 5.0 
# ky   = 0.0
# Epot = 0.0 
# T    = 10.0
# Tout = 0.2


./tdse.x 

# or use the initial data file input

./tdse.x < input 

5. Create animation from generated wave function

./makeanim.sh

6. You can visualize the GIF animate sequence either using display coomand or any web browser

display wf.gif 

or

firefox file://./wf.gif

firefox file://./wf.mp4
