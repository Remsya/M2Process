//
//  StormFunctionLaplace.cpp
//  LaplaceProcess
//
//  Created by Rémi Carnec on 28/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include "StormFunctionLaplace.hpp"

// Constructor
StormFunctionLaplace::StormFunctionLaplace(int dim, double scale, double exponent)
: StormFunction(dim, scale, exponent)
{
    domInfluenceImplemented = (alpha <= 2);

    // Call base class protected method
    computeMoments();
}

// Compute StormFunctionLaplace(x)
double StormFunctionLaplace::computeFunction(double x) const
{
    return exp(-pow(abs(x/scaleFactor),alpha));
}

// Calculate k^th moment of the StormFunctionLaplace
double StormFunctionLaplace::computeMoment(int k) const
{
    return pow(scaleFactor,k+1)*tgamma((k+1)/alpha)/alpha;
}

// TODO : To be converted in enum if really needed
nameShapeFunction StormFunctionLaplace::getShape() const
{
    return nameShapeFunction::Laplace;
}
