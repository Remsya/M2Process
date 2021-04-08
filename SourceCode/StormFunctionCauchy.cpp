//
//  StormFunctionCauchy.cpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 28/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include "StormFunctionCauchy.hpp"

// Constructor
StormFunctionCauchy::StormFunctionCauchy(int dim, double scale, double exponent)
: StormFunction(dim, scale, exponent)
{
    // Check parameters before creating Cauchy Storm Function
    if (exponent < dim/2.)
    {
        throw ("If Storm Function is Cauchy, then alpha must be > dimension/2");
    }
    
    domInfluenceImplemented = true;

    // Call base class protected method
    computeMoments();
}

// Compute StormFunctionCauchy(x)
double StormFunctionCauchy::computeFunction(double x) const
{
    return 1/pow(1+pow(x/scaleFactor, 2),alpha);
}

// Calculate k^th moment of the StormFunctionCauchy
double StormFunctionCauchy::computeMoment(int k) const
{
    return pow(scaleFactor,k+1)*tgamma((k+1.)/2)*tgamma(alpha-(k+1.)/2)/tgamma(alpha)/2;
}

// TODO : To be converted in enum if really needed
nameShapeFunction StormFunctionCauchy::getShape() const
{
    return nameShapeFunction::Cauchy;
}
