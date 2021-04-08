# M2Process

Storm Process simulation algorithm

## Recap

This repository contains the reference implementation of the algorithm detailed in *Continuous simulation of storm processes*, by Demangeot et al., as well as a basic R interface.

## Code
A starting code is provided in ```SimM2Process.R```, where the user can call the C++ source code to simulate within a *d-rectangle* or a *d-sphere*, where *d* can be any dimension. Instructions for 2D-visualization are also provided. For now, two storm processes are available in the library:
- Laplace storms (e.g. Gaussian): ![f(u)=e^{-\left(\dfrac{u}{a}\right)^\alpha}](https://latex.codecogs.com/svg.latex?f(u)=e^{-\left(\dfrac{u}{a}\right)^\alpha})
- Student storms (e.g. Cauchy): ![f(u)=\left(1+\dfrac{u^2}{a^2}\right)^{-\alpha}](https://latex.codecogs.com/svg.latex?f(u)=\left(1+\dfrac{u^2}{a^2}\right)^{-\alpha})

The implementation relies on the following files:
- ```R_SimM2Process.cpp```: Contains the interface between R and C++ to simulate a storm process or compute the maximas when the Poisson points are already calculated and stored.
- ```M2ProcessSimulate.cpp```: Global structure of the simulation of inner and outer processes (see *Continuous simulation of storm processes*).
- ```Grid.cpp```: Precalculates the covering to access the neighbours faster when comparing the domains of influence.
- ```PoissonPoint.cpp```: Poisson point object (with a real and spatial components) which methods allow to compare the domains of influence between two points.
- ```M2ProcessHelper.cpp``` and ```VectorHelper.cpp```: Helper functions to compute necessary mathematical quantities and manipulate vectors.
- ```StormFunction.cpp```, ```StormFunctionLaplace.cpp```, ```StormFunctionCauchy.cpp```: To create and compute the shape functions, as well as their moments.


## Example
The figure below illustrates a simulation of a 2D Student Storm Process in a (400,600) rectangle with ![\alpha=1.5,a=1](https://latex.codecogs.com/svg.latex?\alpha=1.5,a=1).

![alt text](https://github.com/Remsya/M2Process/blob/main/Files/Rplot.png)

