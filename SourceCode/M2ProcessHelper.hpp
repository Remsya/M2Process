//
//  M2ProcessHelper.hpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 21/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//
#ifndef M2ProcessHelper_hpp
#define M2ProcessHelper_hpp


#include <random>

#include "StormFunction.hpp"
#include "Grid.hpp"
#include "VectorHelper.hpp"


// Helper functions to simulate Z1,Z2
double beta(double a, double b);
double comb(int k, int n);
double computeVolumeD1(const StormFunction& shape, Grid *grid, double C_N);
double computeVolumeD2(const StormFunction& shape, Grid *grid, double C_N);
double distanceFromGrid(std::vector<double>& x, Grid *grid);
std::vector<double> weightsAndVolume(const StormFunction& shape, Grid* grid, double& volume, double C);
int countOnes (int n);

// Simulators of random variables
std::vector<double> unifSphere(std::default_random_engine& gen, int dimension, double r);
std::vector<double> unifGrid(std::default_random_engine& gen, Grid *grid);
std::vector<double> simulateDensityBall(std::default_random_engine& gen, const StormFunction& shape, double r);
std::vector<double> simulateDensityRectangle(std::default_random_engine& gen, const StormFunction& shape, Grid *grid, const std::vector<double>& p_I);
double simulateFk(std::default_random_engine& gen, const StormFunction& shape, int k);
double generalizedGammaDist(std::default_random_engine& gen, double a, double beta, double alpha);

#endif /* M2ProcessHelper_hpp */
