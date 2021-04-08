//
//  StormFunction.hpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 28/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#ifndef StormFunctionCauchy_hpp
#define StormFunctionCauchy_hpp

#include "StormFunction.hpp"

// This class represent a cauchy storm function
class StormFunctionCauchy : public StormFunction {
    
public:
    // Constructor
    StormFunctionCauchy(int dimension, double scaleFactor, double alpha);
    
    // Compute StormFunction(x)
    double computeFunction(double x) const override;
    
    // Compute k^th moment
    double computeMoment(int k) const override;

    // Virtual Accessors
    nameShapeFunction getShape() const override;
};

#endif /* StormFunctionCauchy_hpp */
