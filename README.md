# M2Process

Storm Process continuous simulation algorithm.

## Recap

This repository contains the reference implementation of the algorithm detailed in *Continuous simulation of storm processes*, by Demangeot et al., as well as a basic R interface.

## Code
A starting code is provided in ```SimM2Process.R```, where the user can call the C++ source code to simulate within a *d-rectangle* or a *d-sphere*, where *d* can be any dimension (Note that the path should be adapted at the beginning of the file). Instructions for 2D-visualization are also provided. For now, two storm processes are available in the library:
- Laplace storms (e.g. Gaussian): ![f(u)=e^{-\left(u/a\right)^\alpha}](https://latex.codecogs.com/svg.latex?f(u)=e^{-\left(u/a\right)^\alpha})
- Student storms (e.g. Cauchy): ![f(u)=\left(1+u^2/a^2\right)^{-\alpha}](https://latex.codecogs.com/svg.latex?f(u)=\left(1+u^2/a^2\right)^{-\alpha})

The implementation of the algorithm relies on the following files:
- ```R_SimM2Process.cpp```: Contains the interface between R and C++ to simulate a storm process or compute the maximas when the Poisson points are already calculated and stored.
- ```M2ProcessSimulate.cpp```: Global structure of the simulation of the inner and outer processes (see *Continuous simulation of storm processes* for more details).
- ```Grid.cpp```: Precalculates the covering to access the neighbours faster when comparing the domains of influence of two Poisson points.
- ```PoissonPoint.cpp```: Poisson point object (with a real and spatial components) which methods allow to compare the domains of influence between two points if implemented for the corresponding shape function (i.e. to update the box containing the domain of influence of a new point).
- ```M2ProcessHelper.cpp``` and ```VectorHelper.cpp```: Helper functions to compute necessary mathematical quantities (volumes, etc.) and manipulate vectors.
- ```StormFunction.cpp```, ```StormFunctionLaplace.cpp```, ```StormFunctionCauchy.cpp```: To create and compute the shape functions, as well as their moments.

Once the simulation is complete, the Poisson points and their domains of influence are stored in text files (`T.txt`, `S.txt`, `R.txt`) in the folder *Files* so that they can be used to compute the Storm process at any given point within the simulation window (ball or rectangle).

## Example
The figure below illustrates a simulation of a 2D Student Storm Process in a (400,600) rectangle with ![\alpha=1.5,a=1](https://latex.codecogs.com/svg.latex?\alpha=1.5,a=1). Note that this is the default simulation example given in ```SimM2Process.R```, also represented in *Files/CauchySim.pdf*.

![alt text](https://github.com/Remsya/M2Process/blob/main/Files/CauchySim.png)

