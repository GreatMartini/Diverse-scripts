This Jupyter notebook contains a project realised in partnership with my fellow student Alexandre Bougakov. In the project the heat transfer equation is solved numerically using the Crank-Nicholson algorythm. The aim was to determine the heat transfer in a telluric celestial body. 

First the code solves the heat equation in a rod with a sinusoidal function as initial condition. Knowing the anlytical solution, the results are compared to see how the numerical solution presents errors at the borders of the rod and increase with the simulation time.

Afterwards, the code solves the heat equation for a homogeneous spherical body, using the moon's physical properties as parameters.

Then, the code solves the heat equation with a radioactive energy source term while exploring the parameters that define the thermal behaviour of the body. First the physical properties of the moon are used. Afterwards, the code was applyed to the properties of Trans-Npetunian Objects hence, the analysed body was decomposed in one ice layer and a silicate nucleus. The paremeters used were the physical properties of Quaoar.