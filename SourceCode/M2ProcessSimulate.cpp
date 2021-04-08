//
//  M2ProcessSimulate.cpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 21/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include "M2ProcessSimulate.hpp"

// Define simulation functions


bool comparator(PoissonPoint& p1, PoissonPoint& p2)
{
    return p1.getT() < p2.getT();
}

// Simulation in B(0,r)
simulationReturn simulateZ1(std::default_random_engine& gen, const StormFunction& shape, Grid* grid)
{
    // To store poisson points
    simulationReturn simulationResult;
    std::vector<PoissonPoint> z;
    
    // Compute volume of B^d(0,r) and integral of function
    int dim = shape.getDimension();
    double vol = grid->computeVolume();
    double mu = 1/(dim*pow(M_PI, dim/2.)/tgamma(dim/2.+1)*shape.getMoment(dim-1));
    
    // To simulate the poisson points
    std::exponential_distribution<double>expDistribution(1);
    double cumsum = 0;
    std::vector<double> unif;
    
    // To compare domains of influence
    std::vector<std::vector<int>> neighboursIndices = grid->getNeighboursIndices();
    std::vector<std::vector<PoissonPoint>> points(neighboursIndices.size());
    int nbSimulation = 0;
    int nbComparison = 0;
    
    // Fill the covering domains
    int nbDomains = 0;
    while (nbDomains < grid->getNumDomains() && z.size() < 1e8)
    {
        // Generate new Poisson Point
        cumsum += expDistribution(gen);
        unif = unifGrid(gen, grid);
        PoissonPoint p(cumsum/(mu*vol), unif, grid);
        
        // Check if new poisson point has influence
        bool hasInfluence = true;
        if (shape.isDomInfluenceImplemented())
        {
            int index = grid->findIndex(unif);
            for (auto it = neighboursIndices[index].begin(); hasInfluence && it != neighboursIndices[index].end(); ++it)
            {
                int k = 0;
                while (hasInfluence && k < points[*it].size())
                {
                    hasInfluence = p.hasInfluence(points[*it][k], shape);
                    ++k;
                }
                nbComparison += k;
            }
            if (hasInfluence || points[index].empty())
            {
                points[index].push_back(p);
            }
        }
        if (hasInfluence)
        {
            z.push_back(p);
        }
        
        // Check how many balls contain this Poisson Point and did not contain any point before
        nbDomains += grid->findNbDomains(unif);
        ++nbSimulation;
    }
    
    // Save T_N for outer process
    simulationResult.T_N = cumsum/(mu*vol);

    // Simulate until stopping rule occurs
    double upperBound = shape.computeFunction(grid->getDelta())/cumsum;
    while (1./cumsum > upperBound && z.size() < 1e8)
    {
        // Generate new Poisson Point
        cumsum += expDistribution(gen);
        unif = unifGrid(gen, grid);
        PoissonPoint p(cumsum/(mu*vol), unif, grid);

        // Check if new poisson point has influence
        bool hasInfluence = true;
        if (shape.isDomInfluenceImplemented())
        {
            int index = grid->findIndex(unif);
            for (auto it = neighboursIndices[index].begin(); hasInfluence && it != neighboursIndices[index].end(); ++it)
            {
                int k = 0;
                while (hasInfluence && k < points[*it].size())
                {
                    hasInfluence = p.hasInfluence(points[*it][k], shape);
                    ++k;
                }
                nbComparison += k;
            }
            if (hasInfluence)
            {
                points[index].push_back(p);
            }
        }
        if (hasInfluence)
        {
            z.push_back(p);
        }
        
        ++nbSimulation;
    }

    std::cout<<"Inner Process: " << z.size() << " Poisson points insead of: " << nbSimulation << std::endl;
    
    // Store Poisson points
    simulationResult.poissonPoints = z;
    simulationResult.points = points;
    return simulationResult;
}

// Simulation in B^c(0,r)
std::vector<PoissonPoint> simulateZ2(std::default_random_engine& gen, const StormFunction& shape, Grid *grid, double C, std::vector<std::vector<PoissonPoint>>& points)
{
    // To store poisson points
    std::vector<PoissonPoint> z;
    
    // Compute mu
    int dim = shape.getDimension();
    double mu = 1/(dim*pow(M_PI, dim/2.)/tgamma(dim/2.+1)*shape.getMoment(dim-1));
    
    if (grid->getTypeGrid() == TypeGrid::Ball)
    {
        // Distributions
        double vol = computeVolumeD2(shape, grid, C);
        std::poisson_distribution<int>poissonDistribution(mu*vol);
        std::uniform_real_distribution<double>unifDistribution(0,1);

        // Number of points to simulate in D2
        int N = poissonDistribution(gen);
        for (int i = 0; i < N; ++i)
        {
            // Simulate a Poisson point (T,S) uniformly distributed on the desired region
            std::vector<double> S = simulateDensityBall(gen, shape, grid->getRadius());
            double T = C * shape.computeFunction(norm(S)-grid->getRadius()) * unifDistribution(gen);
            
            // Store the corresponding PoissonPoint
            PoissonPoint p(T, S, grid);
            z.push_back(p);
        }
    }
    // If rectangle, simulate in D1 and D2, else, only simulate in D2
    else if (grid->getTypeGrid() == TypeGrid::Rectangle)
    {
        // Build p_I
        double vol = 0;
        std::vector<double> p_I = weightsAndVolume(shape, grid, vol, C);
        
        // Distributions
        std::poisson_distribution<int>poissonDistribution(mu*vol);
        std::uniform_real_distribution<double>unifDistribution(0,1);
        
        // Number of points to simulate in D1
        int N = poissonDistribution(gen);
        for (int i = 0; i < N; ++i)
        {
            // Simulate a point (T,S) uniformly distributed on the desired region
            std::vector<double> S = simulateDensityRectangle(gen, shape, grid, p_I);
            double T = C * unifDistribution(gen) * shape.computeFunction(distanceFromGrid(S, grid));
 
            // Store the corresponding PoissonPoint
            PoissonPoint p(T, S, grid);
            z.push_back(p);
        }
        
    }
    
    // Sort vector and compare to eliminate useless candidates
    if (shape.isDomInfluenceImplemented())
    {
        // Sort Poisson points comparing their T values
        sort(z.begin(), z.end(), &comparator);
        std::vector<PoissonPoint> zNew;
        std::vector<std::vector<int>> neighboursIndices = grid->getNeighboursIndices();
        
        // For each poisson point p, compare it to its neighbours
        for (PoissonPoint& p: z)
        {
            bool hasInfluence = true;
            int index = grid->findNearestDomain(p.getS());
            std::vector<int> indices = neighboursIndices[index];
            for (auto it = indices.begin(); hasInfluence && it != indices.end(); ++it)
            {
                int j = 0;
                while (hasInfluence && j < points[*it].size())
                {
                    hasInfluence = p.hasInfluence(points[*it][j], shape);
                    ++j;
                }
            }
            if (hasInfluence)
            {
                points[index].push_back(p);
                zNew.push_back(p);
            }
        }
        std::cout<<"Outer Process: " << zNew.size() << " Poisson points instead of: " << z.size() << std::endl;
        return zNew;
    }
    else
    {
        std::cout<<"Outer Process: " << z.size() << " Poisson points." << std::endl;
        return z;
    }
    
}
