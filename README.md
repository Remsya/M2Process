# M2Process

Storm Process simulation algorithm

## Recap

This repository contains the reference implementation of the algorithm detailed in *Continuous simulation of storm processes*, by Demangeot et al., as well as a basic R interface.

## Code
A starting code is provided in ```SimM2Process.R```, where the user can call the C++ source code to simulate within a *d-rectangle* or a *d-sphere*, where *d* can be any dimension. For *d = 2*, a visualization code is provided. For now, two storm processes are available:
- Laplace storms: ![x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}]
- Cauchy

![alt text](https://github.com/Remsya/M2Process/blob/main/Files/Rplot.png)

A starting code is provided in ```main_py```, where the user can play around with the different methods on two datasets that you will find in the `data` folder: the *Stanford Bunny* and a lighter version of the original *Asian Dragon*, downloaded from the Standord 3D Scanning Repository (https://graphics.stanford.edu/data/3Dscanrep/). The implementation relies on three main files:
- ```optimize.py```: contains the main structure of the optimization algorithm (see *Generalized-ICP*, Segal et al.).
- ```algorithms.py```: used to minimize the loss at each iteration for the different methods (using a closed form solution for *Standard ICP* and the Conjugate Gradient for *Point-to-plane* and *Generalized-ICP*).
- ```utils.py```: contains all sorts of useful functions to perform PCA, compute rotation matrices, gradients, etc.
