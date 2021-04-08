#### Example of simulation of the M2 field
#### Vizualization

#### Author: Rémi Carnec
#### Date: 06/07/2020

### Link with libraries

# Imports
rm(list = ls())
library(Rcpp)
library(RcppArmadillo)

# Link with C++ files. IMPORTANT: NEED TO SET 'path' VARIABLE BEFORE RUNNING SIMULATION
path <- "Desktop/"
SaveFilesPath <- paste(path, "RM2Process/Files/", sep = "")
sourceCpp(file = paste(path, "RM2Process/SourceCode/R_SimM2Process.cpp", sep=""))


### Simulation M2 Process

# Storm function params
shape <- "Cauchy"
alpha <- 3./2
a <- 1
dim <- 2

# Window and covering params (for rectangle)
window <- c(100,200)
covering <- rep(1, dim)

# Window and covering params (for ball)
radius <- 20
delta <- 1

# Build grid of points (if ball, grid must be included in B(0,radius)). Does not have to be regular.
nPoints <- c(100,200)
Grid <- generateNRectangle(window, nPoints)

# Simulate M2 Process. If seed is not defined, then the process is random. simulation[i] is the value in Grid[i].
simulation <- R_M2ProcessRectangle(shape, dim, alpha, a, window, covering, Grid, path = SaveFilesPath, seed = 1)

# Convert Data to matrix to plot
simulation_mat <- matrix(0, nrow = nPoints[1], ncol = nPoints[2])
for (i in 1:nPoints[1])
{
  for (j in 1:nPoints[2])
  {
    simulation_mat[i,j] <- simulation[i+(j-1)*nPoints[1]]
  }
}

# Color package
#install.packages("viridis")
#library("viridis")

# Convert grid to plot. Need the grid to be regular.
v_x <- seq(-window[1]/2, window[1]/2, length.out = nPoints[1])
v_y <- seq(-window[2]/2, window[2]/2, length.out = nPoints[2])

# Vizualization of log(Max stable function). ONLY FOR 2D.
pdf(file = paste(SaveFilesPath,"Maxstab.pdf", sep =""), width = 6.8, height = 6)
filled.contour(v_x, v_y, log(simulation_mat), color.palette = viridis)
dev.off()


### If need to evaluate the function in other points using Poisson points values

# Get Poisson points data of previous simulation (T, S, rect of influence)
T_vect <- as.vector(read.delim(paste(SaveFilesPath,"T.txt", sep =""), header = F)[[1]])
S_raw <-as.vector(read.delim(paste(SaveFilesPath,"S.txt", sep =""), header = F)[[1]])
Rect_raw <- as.vector(read.delim(paste(SaveFilesPath,"R.txt", sep =""), header = F)[[1]])
numPoissonPoints <- length(T_vect)
S_vect <- list()
Rect_vect <- list()
for (i in 1:numPoissonPoints){
  # Build S_i
  S_vect[i] <- list(S_raw[(1+(i-1)*dim):(i*dim)])
  
  # Build Rect_i
  Rect_vect[(i-1)*2+1] <- list(Rect_raw[(2*(i-1)*dim +1) : ((2*i-1)*dim)])
  Rect_vect[i*2] <- list(Rect_raw[((2*i-1)*dim+1):(2*i*dim)])
}

# Evaluate M2 Process in desired points
simulation <- R_EvaluateM2Process(shape, dim, alpha, a, window, covering, Grid, T_vect, S_vect, Rect_vect)
