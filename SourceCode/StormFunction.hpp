//
//  StormFunction.hpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 28/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#ifndef StormFunction_hpp
#define StormFunction_hpp

#include <stdio.h>
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>
#include <iostream>

// This class represent a storm function f : R^d -> R^+, with name "shape", d = dimension

enum nameShapeFunction
{
    Cauchy,
    Laplace
};

class StormFunction {
    
protected:
    int dimension;
    double scaleFactor;
    double alpha;
    double cumulativeMoments;
    std::vector<double> moments;
    bool domInfluenceImplemented;
    
    void computeMoments();
    
public:
    // Constructor
    StormFunction(int dimension, double scaleFactor, double alpha);
    
    // Compute StormFunction(x)
    virtual double computeFunction(double x) const = 0;
    
    // Compute k^th moment
    virtual double computeMoment(int k) const = 0;
    
    // Virtual Accessors
    virtual nameShapeFunction getShape() const = 0;
    
    // Accessors
    inline int getDimension() const
    {
        return dimension;
    }
    
    inline double getMoment(int i) const
    {
        // Check that 0 <= i < dimension
        if (i < 0 || i >= dimension)
        {
            throw("Moment index must be between 0 and dim - 1");
        }
        return moments[i];
    }
    
    inline double getCumMoment() const
    {
        return cumulativeMoments;
    }
    
    inline double getScaleFactor() const
    {
        return scaleFactor;
    }
    
    inline double getAlpha() const
    {
        return alpha;
    }
    
    inline bool isDomInfluenceImplemented() const
    {
        return domInfluenceImplemented;
    }
};

#endif /* StormFunction_hpp */
