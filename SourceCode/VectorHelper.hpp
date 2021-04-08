//
//  VectorHelper.hpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 16/06/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#ifndef VectorHelper_hpp
#define VectorHelper_hpp

#include <stdio.h>
#include <vector>
#include <math.h>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>

// Helper functions to process vectors and data
double norm(const std::vector<double>& x);
std::vector<double> addConstant(const std::vector<double> x, double m);
std::vector<double> linearCombination(const std::vector<double>& x, const std::vector<double>& y, double m_x, double m_y);
std::vector<double> multiplyByScalar(const std::vector<double> x, double m);
std::vector<double> addVector(const std::vector<double>& x, const std::vector<double>& y);
std::vector<double> substractVector(const std::vector<double>& x, const std::vector<double>& y);
double distance(const std::vector<double>& x, const std::vector<double>& y);

#endif /* VectorHelper_hpp */
