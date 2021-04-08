# M2Process

Storm Process simulation algorithm

## Recap

This repository contains the reference implementation of the algorithm detailed in *Continuous simulation of storm processes*, by Demangeot et al., as well as a basic R interface.

## Code
A starting code is provided in ```SimM2Process.R```, where the user can call the C++ source code to simulate within a *d-rectangle* or a *d-sphere*, where *d* can be any dimension. Instructions for 2D-visualization are also provided. For now, two storm processes are available in the library:
- Laplace storms
![f(x)=e^{-\frac{x^2}{a^2}}](https://latex.codecogs.com/svg.latex?f(x)=e^{-\frac{x^2}{a^2}})
- Student storms: 
![f(x)=\frac{1}{\left(1+x^2/a^2\right)^\alpha}](https://latex.codecogs.com/svg.latex?f(x)=\frac{1}{\left(1+x^2/a^2\right)})
![f(u)=\left(1+\dfrac{u^2}{a^2}\right)^{-\alpha}](https://latex.codecogs.com/svg.latex?f(u)=\left(1+\dfrac{u^2}{a^2}\right)^{-\alpha})




A starting code is provided in ```main_py```, where the user can play around with the different methods on two datasets that you will find in the `data` folder: the *Stanford Bunny* and a lighter version of the original *Asian Dragon*, downloaded from the Standord 3D Scanning Repository (https://graphics.stanford.edu/data/3Dscanrep/). The implementation relies on three main files:
- ```optimize.py```: contains the main structure of the optimization algorithm (see *Generalized-ICP*, Segal et al.).
- ```algorithms.py```: used to minimize the loss at each iteration for the different methods (using a closed form solution for *Standard ICP* and the Conjugate Gradient for *Point-to-plane* and *Generalized-ICP*).
- ```utils.py```: contains all sorts of useful functions to perform PCA, compute rotation matrices, gradients, etc.

## Example
The figure below illustrates a simulation of a 2D Student Storm Process in a (400,600) rectangle with ![\alpha=1.5,a=1](https://latex.codecogs.com/svg.latex?\alpha=1.5,a=1).

![alt text](https://github.com/Remsya/M2Process/blob/main/Files/Rplot.png)

