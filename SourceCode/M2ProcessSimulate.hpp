//
//  M2ProcessSimulate.hpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 21/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#ifndef M2ProcessSimulate_hpp
#define M2ProcessSimulate_hpp

#include "PoissonPoint.hpp"
#include "M2ProcessHelper.hpp"
#include "Grid.hpp"
#include "StormFunctionLaplace.hpp"
#include "StormFunctionCauchy.hpp"

struct simulationReturn
{
    // Vector of Poisson Points after simulation of Z1
    std::vector<PoissonPoint> poissonPoints;
    std::vector<std::vector<PoissonPoint>> points;
    
    // Value of T_k when all balls contain at least one Poisson point
    double T_N;
};

simulationReturn simulateZ1(std::default_random_engine& gen, const StormFunction& shape, Grid* grid);
std::vector<PoissonPoint> simulateZ2(std::default_random_engine& gen, const StormFunction& shape, Grid* grid, double C, std::vector<std::vector<PoissonPoint>>& points);

// Main function that split simulation between B(0,r) and B^c(0,r)
template <typename T> std::vector<PoissonPoint> M2Process(std::string shape, int dim, double alpha, double a, T paramGrid, T paramCovering, std::default_random_engine& generator){
    // Generate Grid
    Grid* grid(new Grid(paramGrid, paramCovering, dim));
    
    // Create StormFunction
    StormFunction* shapeFunction;
    if (shape == "Laplace")
    {
        shapeFunction = new StormFunctionLaplace(dim, a, alpha);
    }
    else if (shape == "Cauchy")
    {
        shapeFunction = new StormFunctionCauchy(dim, a, alpha);
    }
    else
    {
        throw("Unknown Storm Function shape");
    }
    
    // Simulating samples inside Grid, returns a first vector of Poisson Points
    simulationReturn simulationResult = simulateZ1(generator, *shapeFunction, grid);
    
    // Calculate C_N parameter for Z2
    double C_N = simulationResult.T_N/shapeFunction->computeFunction(grid->getDelta());
    
    // Simulating samples outside Grid
    std::vector<PoissonPoint> z2 = simulateZ2(generator, *shapeFunction, grid, C_N, simulationResult.points);
    
    // Gather data
    simulationResult.poissonPoints.insert(simulationResult.poissonPoints.end(), z2.begin(), z2.end());
    simulationResult.points.clear();
    
    return simulationResult.poissonPoints;
};

#endif /* M2ProcessSimulate_hpp */
