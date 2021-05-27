# micropipette_aggregate

This program is developed to simulate 2d cell membranes and their interactions.
Membranes are composed by several particles connected.

The membrane is ruled by 4 interactions:
 - Springs between the composing particles;
 - Bending potential to minimise the angle between the particles;
 - Area conservation potential;
 - Linear contact repulsion.
 
The interaction between two distinct cells is given by:
 - Linear contact repulsion;
 - Linear short-range attraction.
 
How the code works:
 - This module should be imported by python, and on its first run should be send the "True" parameter to the init function, to initialize the compiling routine.
 - The time evolution is performed by Fortran95 routine, maximising the code performance.
 - Python send the current cells state to Fortran95 routine, where the particles evolve for several time steps and return the a new particles configuration.
 - The integer parameter given to the integrate Fortran95 function defines the number of cores that would be used in the calculations.
 - For the current time, due to the sequencial regions of the Fortran routine, I strongly recommend to not use more than 3 cores
