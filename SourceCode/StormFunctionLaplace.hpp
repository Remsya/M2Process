//
//  StormFunction.hpp
//  LaplaceProcess
//
//  Created by Rémi Carnec on 28/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#ifndef StormFunctionLaplace_hpp
#define StormFunctionLaplace_hpp

#include "StormFunction.hpp"

// This class represent a laplace storm function
class StormFunctionLaplace : public StormFunction {
    
public:
    // Constructor
    StormFunctionLaplace(int dimension, double scaleFactor, double alpha);
    
    // Compute StormFunction(x)
    double computeFunction(double x) const override;
    
    // Compute k^th moment
    double computeMoment(int k) const override;

    // Virtual Accessors
    nameShapeFunction getShape() const override;
};

#endif /* StormFunctionLaplace_hpp */
