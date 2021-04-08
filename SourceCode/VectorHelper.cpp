//
//  VectorHelper.cpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 16/06/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include "VectorHelper.hpp"

// Helper functions to deal with data and vectors

// Norm of a vector
double norm(const std::vector<double>& x)
{
    double s = 0;
    for (double i: x)
    {
        s += i*i;
    }
    return sqrt(s);
}

std::vector<double> addConstant(std::vector<double> x, double m)
{
    // add constant
    for (double& element: x)
    {
        element += m;
    }
    
    return x;
}

std::vector<double> multiplyByScalar(std::vector<double> x, double m)
{
    std::transform(x.begin(), x.end(), x.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, m));
    return x;
}


std::vector<double> addVector(const std::vector<double>& x, const std::vector<double>& y)
{
    // Check same size
    if (x.size() != y.size())
    {
        throw("Size of x and y incoherent in addVector");
    }
    
    // Build result
    std::vector<double> result = x;
    std::transform(result.begin(), result.end(), y.begin(), result.begin(), std::plus<double>());
    
    return result;
}

std::vector<double> substractVector(const std::vector<double>& x,const std::vector<double>& y)
{
    // Check same size
    if (x.size() != y.size())
    {
        throw("Size of x and y incoherent in substractVector");
    }
    
    // Build result
    std::vector<double> result = x;
    std::transform(result.begin(), result.end(), y.begin(), result.begin(), std::minus<double>());
    
    return result;
}

double distance(const std::vector<double>& x, const std::vector<double>& y)
{
    double dist = 0;
    
    double s;
    for (int i = 0; i < x.size(); ++i)
    {
        s = x[i]-y[i];
        dist += s*s;
    }
    
    return sqrt(dist);
}

std::vector<double> linearCombination(const std::vector<double>& x, const std::vector<double>& y, double m_x, double m_y)
{
    // Check same size
    if (x.size() != y.size())
    {
        throw("Size of x and y incoherent in linearCombination");
    }
    
    std::vector<double> result(x.size(), 0.);
    for (int i = 0; i < x.size(); ++i)
    {
        result[i] = m_x * x[i] + m_y * y[i];
    }
    return result;
}
