//
//  StormFunction.cpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 28/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include "StormFunction.hpp"

// Constructor
StormFunction::StormFunction(int dim, double scale, double exponent)
{
    // Check parameters before creating Storm Function
    if (dim < 1)
    {
        throw("Dimension must be at least 1");
    }
    if (scale <= 0)
    {
        throw("Scaling factor must be > 0");
    }
    
    // Initialize StormFunction parameters
    dimension = dim;
    scaleFactor = scale;
    alpha = exponent;
    domInfluenceImplemented = false;
}

// Calculate k^th moment of the StormFunction
void StormFunction::computeMoments()
{
    // Compute moments to avoid recalculating them
    cumulativeMoments = 0;
    moments.resize(dimension);
    for (int i = 0; i < dimension; ++i)
    {
        moments[i] = computeMoment(i);
        cumulativeMoments += moments[i];
    }
}
