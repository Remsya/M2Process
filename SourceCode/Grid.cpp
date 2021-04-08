//
//  Grid.cpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 03/06/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include "Grid.hpp"

// Constructor for grid in Ball
Grid::Grid(double R, double Delta, int Dimension)
{
    // Assign values
    r = R;
    delta = Delta;
    dimension = Dimension;
    typeGrid = TypeGrid::Ball;
    
    // Number of centers per axis
    nAxis = 2*int(2*r/delta)+1;
    
    // Initialize activation vector: activatedBalls[i] = true if i_th ball is in B(0,r) and does not contain a Poisson Point yet
    activatedDomains = std::vector<bool>(pow(nAxis, dimension), false);
    
    // Number of points contained in the ball
    nbDomains = 0;
    
    // Build coordinates of the centers of balls
    grid = {};
    for (int k = 0; k < pow(nAxis,dimension); ++k)
    {
        grid.push_back({});
        for (int i = 0; i < dimension; ++i)
        {
            grid[k].push_back(delta/2.*(int(k*1./pow(nAxis,i))%nAxis-(nAxis-1)/2.));
        }
        
        // If x[k] is in the ball, increment nBall
        if (norm(grid[k])<=r)
        {
            activatedDomains[k] = true;
            ++nbDomains;
        }
    }
    
    // Build neighbours indices
    neighboursIndices.resize(pow(nAxis, dimension));
    int n = 2;
    int nTot = 2*n+1;
    for (int k = 0; k < pow(nAxis, Dimension); ++k)
    {
        // Compute indices of nearest neigbours
        if (activatedDomains[k])
        {
            int nNeighbours = pow(nTot, Dimension);
            for (int i = 0; i < nNeighbours; ++i)
            {
                bool isInGrid = true;
                int index = k;
                for (int j = 0; isInGrid && j < Dimension; ++j)
                {
                    int index_j = int(i/pow(nTot,j))%nTot - n;
                    isInGrid = isInGrid && (index_j == 0 || ((index_j + int(k/pow(nAxis,j))%nAxis >= 0)
                            && (int(k/pow(nAxis,j))%nAxis + index_j <= nAxis-1)));
                    index += index_j*pow(nAxis,j);
                }
                if (isInGrid && norm(grid[index]) <= r)
                {
                    neighboursIndices[k].push_back(index);
                }
            }
        }
    }
}

// Constructor for grid in Rectangle
Grid::Grid(std::vector<double> SidesRect, std::vector<double> SidesDomains, int Dimension)
{
    // Check sizes are coherent with dimension
    if (SidesDomains.size() != Dimension || SidesRect.size() != Dimension)
    {
        throw("Dimensions of window of simulation or covering different from dimension of Space");
    }

    // Assign values
    typeGrid = TypeGrid::Rectangle;
    dimension = Dimension;
    sidesRect = SidesRect;
    sidesDomains = SidesDomains;
    r = norm(sidesRect)/2;
    delta = norm(sidesDomains);

    // Compute number of rectangles
    nbDomains = round(std::accumulate(sidesRect.begin(), sidesRect.end(), 1., std::multiplies<double>())/std::accumulate(sidesDomains.begin(), sidesDomains.end(), 1., std::multiplies<double>()));
    
    // Initialize activation vector: activatedBalls[i] = true if the i_th domain does not contain a Poisson Point yet
    activatedDomains = std::vector<bool>(nbDomains, true);
    
    // Build coordinates of the centers of the rectangles
    neighboursIndices.resize(nbDomains);
    grid = {};
    
    // Cumproduct
    std::vector<double> cumprod(dimension+1, 1.);
    for (int i = 0; i < dimension; ++i)
    {
        cumprod[i+1] = cumprod[i]*sidesRect[i]/sidesDomains[i];
    }
    
    
    for (int k = 0; k < nbDomains; ++k)
    {
        grid.push_back({});
        for (int i = 0; i < dimension; ++i)
        {
            // Center of rectangle
            int k_i = int(k/cumprod[i])%int(sidesRect[i]/sidesDomains[i]);
            grid[k].push_back((sidesDomains[i]-sidesRect[i])/2+k_i*sidesDomains[i]);
        }
        
        
        // Compute indices of nearest neigbours
        int n = 2;
        int nTot = 2*n+1;
        int nNeighbours = pow(nTot, Dimension);
        for (int i = 0; i < nNeighbours; ++i)
        {
            bool isInGrid = true;
            int index = k;
            for (int j = 0; isInGrid && j < Dimension; ++j)
            {
                int index_j = int(i/pow(nTot,j))%nTot - n;
                isInGrid = isInGrid && (index_j == 0 || ((index_j + int(k/cumprod[j])%int(sidesRect[j]/sidesDomains[j]) >= 0)
                        && (int(k/cumprod[j])%int(sidesRect[j]/sidesDomains[j]) + index_j <= int(sidesRect[j]/sidesDomains[j])-1)));
                index += index_j*cumprod[j];
            }
            if (isInGrid)
            {
                neighboursIndices[k].push_back(index);
            }
        }
        
    }
}

// Accessors
TypeGrid Grid::getTypeGrid()
{
    return typeGrid;
}

double Grid::getRadius()
{
    return r;
}

double Grid::getDelta()
{
    return delta;
}

int Grid::getDimension()
{
    return dimension;
}

int Grid::getNumDomains()
{
    return nbDomains;
}

std::vector<std::vector<double>> Grid::getGrid()
{
    return grid;
}

std::vector<double> Grid::getSidesRect()
{
    return sidesRect;
}

std::vector<double> Grid::getSidesDomains()
{
    return sidesDomains;
}

std::vector<bool> Grid::getActivatedDomains()
{
    return activatedDomains;
}

std::vector<std::vector<int>> Grid::getNeighboursIndices()
{
    return neighboursIndices;
}

// Compute the volume of the grid (either ball or rectangle)
double Grid::computeVolume()
{
    double volume = 0;
    
    if (typeGrid == TypeGrid::Ball)
    {
        volume = pow(r, dimension)*pow(M_PI, dimension/2.)/tgamma(dimension/2.+1);
    }
    else if (typeGrid == TypeGrid::Rectangle)
    {
        volume = std::accumulate(sidesRect.begin(), sidesRect.end(), 1, std::multiplies<double>());
    }
    
    return volume;
}

// Method to find the index of a given point of the grid
int Grid::findIndex(std::vector<double> x)
{
    double index = 0;
    if (typeGrid == TypeGrid::Ball)
    {
        for (int i = 0; i < dimension; ++i)
        {
            index += pow(nAxis,i) * round(2*x[i]/delta+(nAxis-1)/2.);
        }
    }
    else if (typeGrid == TypeGrid::Rectangle)
    {
        double cumProd = 1;
        for (int i = 0; i < dimension; ++i)
        {
            if (abs(x[i]) == sidesRect[i]/2)
            {
                x[i] += (x[i] > 0? -1 : 1) * sidesDomains[i]/2;
            }
            index += int((x[i]+sidesRect[i]/2)/sidesDomains[i])*cumProd;
            cumProd *= sidesRect[i]/sidesDomains[i];
        }
    }
    
    return index;
}

// Method to find all indices of domains that contain a point x
std::vector<int> Grid::findIndices(const std::vector<double>& x)
{
    std::vector<int> indices;
    if (typeGrid == TypeGrid::Ball)
    {
        int nDomains = pow(2,dimension);
        for (int i = 0; i < nDomains; ++i)
        {
            int index = 0;
            for (int j = 0; j < dimension; ++j)
            {
                int indicator_j = int(i/pow(2.,j))%2;
                index += pow(nAxis,j)*int(2*x[j]/delta+indicator_j+(nAxis-1)/2.);
            }
            // Problem of accuracy in C++, need to round the index
            indices.push_back(round(index));
        }
        
    }
    else if (typeGrid == TypeGrid::Rectangle)
    {
        double cumProd = 1;
        int index = 0;
        for (int i = 0; i < dimension; ++i)
        {
            index += int((x[i]+sidesRect[i]/2)/sidesDomains[i])*cumProd;
            cumProd *= sidesRect[i]/sidesDomains[i];
        }
        indices.push_back(index);
    }
    
    return indices;
}

// Method to find how many balls contain x, and update activatedBalls
int Grid::findNbDomains(const std::vector<double>& x)
{
    // Domain Counter
    int n = 0;
    std::vector<int> indices = findIndices(x);
    
    if (typeGrid == TypeGrid::Ball)
    {
        for (int index: indices)
        {
           // If x is in ith Ball AND is the first point to be, increase counter
            if (activatedDomains[index] == true && norm(substractVector(x, grid[index])) <= delta/2.)
            {
                activatedDomains[index] = false;
                ++n;
            }
        }
    }
    else if (typeGrid == TypeGrid::Rectangle)
    {
        for (int index: indices)
        {
            if (activatedDomains[index] == true)
            {
                activatedDomains[index] = false;
                ++n;
            }
        }
    }
    
    return n;
}

// Find covering domain that is closest to x
int Grid::findNearestDomain(const std::vector<double>& x)
{
    std::vector<double> proj(dimension,0.);
    
    if (typeGrid == TypeGrid::Ball)
    {
        // Project x on sphere of radius r
        proj = multiplyByScalar(x, r/norm(x));
    }
    else if (typeGrid == TypeGrid::Rectangle)
    {
        // Rectangle to project
        std::vector<double> lengths = sidesRect;
        
        for (int i = 0; i < dimension; ++i)
        {
            if (x[i] > lengths[i]/2)
            {
                proj[i] = lengths[i]/2-1./(1e9);
            }
            else if (x[i] < -lengths[i]/2)
            {
                proj[i] = -lengths[i]/2+1./(1e9);
            }
        }
    }
    
    return findIndex(proj);
}
