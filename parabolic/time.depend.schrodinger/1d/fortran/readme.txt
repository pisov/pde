# This example demonstrate time evolution of wave packet
#
# Compile the source code:

gfortran -O3 tdse.f90 -o tdse.x -llapack

# Clear previous data files if any

rm -f wf*

# Execute the program to evolve the wave function in time
# Exeample input values
# dx   = 0.01
# dt   = 0.001
# Ekin = 5.0 
# Epot = 10.0
# T    = 8.0
# Tout = 0.1

./tdse.x > energy.dat 

# or use the initial data file input

./tdse.x < input > energy.dat

# Create animation from generated wave function

./makeanim.sh

# You can visualize the GIF animate sequence either using display coomand or any web browser

display wf.gif 

or

firefox file://./wf.gif

firefox file://./wf.mp4
