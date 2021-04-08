//
//  R_SimM2Process.cpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 10/07/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include <R.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <chrono>
#include <fstream>
#include <iterator>
#include <thread>
#include "M2ProcessSimulate.hpp"
#include "M2ProcessHelper.hpp"

using namespace std::chrono;

// If use serial buildMaxFunction, can display estimated time
bool displayTime = false;

// Random number generator
std::random_device rd;

// Files to save data
std::string TFile = "T.txt";
std::string SFile = "S.txt";
std::string RFile = "R.txt";

// Function to save Poisson Points
void savePoissonPoints(std::vector<PoissonPoint>& z, std::string path)
{
    if (z.size() != 0)
    {
        // Output files
        std::ofstream output_T(path + TFile);
        std::ofstream output_S(path + SFile);
        std::ofstream output_Rect(path + RFile);

        // Set precision
        output_T.precision(12);
        output_S.precision(12);
        output_Rect.precision(12);
        
        // Vector iterator
        std::ostream_iterator<double> output_S_it(output_S, "\n");
        std::ostream_iterator<double> output_Influence_it(output_Rect, "\n");
        
        // Save all poisson points
        for (PoissonPoint& p: z)
        {
            // Save T
            output_T << p.getT() << std::endl;
            
            // Save S
            std::copy(p.getS().begin(), p.getS().end(), output_S_it);
            
            // Save domain of influence
            const std::vector<double>& corner1 = p.getRectInfluence()[0];
            const std::vector<double>& corner2 = p.getRectInfluence()[1];
            std::copy(corner1.begin(), corner1.end(), output_Influence_it);
            std::copy(corner2.begin(), corner2.end(), output_Influence_it);
        }
        
        // Close files
        output_T.close();
        output_S.close();
        output_Rect.close();
    }
}

// Generate grid
// [[Rcpp::export]]
std::vector<std::vector<double>> generateNRectangle(std::vector<double> sidesRect, std::vector<int> nPoints)
{
    // Dimension
    int dimension = sidesRect.size();
    
    // Nb points
    int nbPoints = std::accumulate(nPoints.begin(), nPoints.end(), 1, std::multiplies<int>());
    
    // Compute cumulative products
    std::vector<int> prodNum(dimension+1, 1);
    for (int i = 1; i < prodNum.size(); ++i)
    {
        prodNum[i] = prodNum[i-1]*nPoints[i-1];
    }

    // Compute grid
    std::vector<std::vector<double>> grid(nbPoints);
    for (int k = 0; k < nbPoints; ++k)
    {
        for (int i = 0; i < dimension; ++i)
        {
            // Center of rectangle
            int k_i = int(k/prodNum[i])%nPoints[i];
            grid[k].push_back(-sidesRect[i]/2 + k_i*sidesRect[i]*1./(nPoints[i]-1));
        }
    }
    return grid;
}

// Build max function - Serial version
std::vector<double> buildMaxFunction(const std::vector<PoissonPoint>& poissonPoints, const StormFunction& f, const std::vector<std::vector<double>>& x)
{
    // Store max function
    int xsize = x.size();
    std::vector<double> Y(xsize);
    
    // Display execution time
    if (displayTime)
    {
        std::cout << "Build Max function progress: "<< std::endl;
    }
    auto start = high_resolution_clock::now();
    int i = 0;
    
    // Find max values for each x_i
    for (const std::vector<double>& xpt : x)
    {
        double temp = 0;
        for (const PoissonPoint& p: poissonPoints)
        {
            if (p.contributes(xpt))
            {
                temp = std::max(temp, f.computeFunction(distance(p.getS(), xpt))/p.getT());
            }
        }
        if (displayTime && i % 1000 == 0)
        {
            // Compute time left
            auto duration = duration_cast<milliseconds>(high_resolution_clock::now() - start);
            double time_left = (i == 0) ? duration.count()*xsize : duration.count()*(xsize-i)/i;
            std::cout << 100.*(i+1.)/xsize << "%, estimated time left: " << int(time_left/1000) << " s" << std::endl;
        }
        Y[i] = temp;
        ++i;
    }
    return Y;
}

// Build max function - Parallel version
std::vector<double> parallelBuildMaxFunction(const std::vector<PoissonPoint>& poissonPoints, const StormFunction& f, const std::vector<std::vector<double>>& x)
{
    // Store max function
    std::size_t xsize = x.size();
    std::vector<double> Y(xsize);

    // Function to compute max on part of Grid
    auto computeMax = [] (const std::vector<PoissonPoint>& poissonPointsArg, const std::vector<std::vector<double>>& xArg, const StormFunction& fArg, std::vector<double>& YArg, int startIndex, int endIndex)
    {
        // loop over all elements of subvector
        for(int i = startIndex; i < endIndex; ++i)
        {
            double temp = 0;
            for (const PoissonPoint& p: poissonPointsArg)
            {
                if (p.contributes(xArg[i]))
                {
                    temp = std::max(temp, fArg.computeFunction(distance(p.getS(), xArg[i]))/p.getT());
                }
            }
            YArg[i] = temp;
        }
    };
    
    // number of threads
    const size_t nthreads = std::thread::hardware_concurrency();
    std::cout<<"Parallel ("<<nthreads<<" threads):"<<std::endl;
    
    // Create threads
    std::vector<std::thread> threads(nthreads);
    for(int t = 0; t < nthreads; t++)
    {
        threads[t] = std::thread(std::ref(computeMax), std::ref(poissonPoints), std::ref(x), std::ref(f), std::ref(Y), t*xsize/nthreads, (t+1) == nthreads? xsize : (t+1)*xsize/nthreads);
    }
    
    // Iterate on threads
    std::for_each(threads.begin(), threads.end(), [](std::thread& t)
    {
        t.join();
    });
    
    return Y;
}

// [[Rcpp::export]]
std::vector<double> R_M2ProcessBall(std::string shape, int dim, double alpha, double a, double paramGrid, double paramCovering, std::vector<std::vector<double>> grid, std::string path, int seed = -1)
{
    // Random number generator
    seed = seed == -1? rd() : seed;
    std::default_random_engine generator(seed);
    
    // Build corresponding shape function
    StormFunction* f;
    if (shape == "Laplace")
        f = new StormFunctionLaplace(dim, a, alpha);
    else if (shape == "Cauchy")
        f = new StormFunctionCauchy(dim, a, alpha);
    else
        throw("Unknown Storm Function shape!");
    
    // Simulate Poisson points
    auto start = high_resolution_clock::now();
    std::vector<PoissonPoint> z = M2Process<double>(shape, dim, alpha, a, paramGrid, paramCovering, generator);
    std::cout << "Simulation of Poisson points took: " << duration_cast<milliseconds>(high_resolution_clock::now() - start).count()
        << " ms" << std::endl;
        
    // Evaluate function in Grid
    start = high_resolution_clock::now();
    std::vector<double> Y = parallelBuildMaxFunction(z, *f, grid);
    std::cout << "Computation of maximas took: " << duration_cast<milliseconds>(high_resolution_clock::now() - start).count()
    << " ms" << std::endl;
    
    // Save Poisson points
    savePoissonPoints(z, path);
    
    return Y;
}

// [[Rcpp::export]]
std::vector<double> R_M2ProcessRectangle(std::string shape, int dim, double alpha, double a, std::vector<double> paramGrid, std::vector<double> paramCovering, std::vector<std::vector<double>> grid, std::string path, int seed = -1)
{
    // Random number generator
    seed = seed == -1? rd() : seed;
    std::default_random_engine generator(seed);
    
    // Build corresponding shape function
    StormFunction* f;
    if (shape == "Laplace")
        f = new StormFunctionLaplace(dim, a, alpha);
    else if (shape == "Cauchy")
        f = new StormFunctionCauchy(dim, a, alpha);
    else
        throw("Unknown Storm Function shape!");
    
    // Simulate Poisson points
    auto start = high_resolution_clock::now();
    std::vector<PoissonPoint> z = M2Process<std::vector<double>>(shape, dim, alpha, a, paramGrid, paramCovering, generator);
    std::cout << "Simulation of Poisson points took: " << duration_cast<milliseconds>(high_resolution_clock::now() - start).count()
        << " ms" << std::endl;
        
    // Evaluate function in Grid
    start = high_resolution_clock::now();
    std::vector<double> Y = parallelBuildMaxFunction(z, *f, grid);
    std::cout << "Computation of maximas took: " << duration_cast<milliseconds>(high_resolution_clock::now() - start).count()
    << " ms" << std::endl;
    
    // Save Poisson points
    savePoissonPoints(z, path);
    
    return Y;
}

// [[Rcpp::export]]
std::vector<double> R_EvaluateM2Process(std::string shape, int dim, double alpha, double a, std::vector<double> paramGrid, std::vector<double> paramCovering, std::vector<std::vector<double>> grid, std::vector<double> T, std::vector<std::vector<double>> S, std::vector<std::vector<double>> rectInfluence)
{
    // Some checks
    if (S.empty() || S.size() != T.size() || 2*S.size() != rectInfluence.size())
    {
        throw("S, T or rectInfluence does not have the right size");
    }
    
    // Create Poisson points
    std::vector<PoissonPoint> z(S.size());
    for (int i = 0; i < S.size(); ++i)
    {
        // Check
        if (S[i].size() != dim || rectInfluence[2*i].size() != dim || rectInfluence[2*i+1].size() != dim)
        {
            throw("S or rectInfluence vectors do not have the right dimension");
        }
        
        // Create Poisson Point
        PoissonPoint p(T[i], S[i], rectInfluence[2*i], rectInfluence[2*i+1]);
        z[i] = p;
    }
    
    // Build corresponding shape function
    StormFunction* f;
    if (shape == "Laplace")
        f = new StormFunctionLaplace(dim, a, alpha);
    else if (shape == "Cauchy")
        f = new StormFunctionCauchy(dim, a, alpha);
    else
        throw("Unknown Storm Function shape!");
    
    // Evaluate function in Grid
    auto start = high_resolution_clock::now();
    std::vector<double> Y = parallelBuildMaxFunction(z, *f, grid);
    std::cout << "Computation of maximas took: " << duration_cast<milliseconds>(high_resolution_clock::now() - start).count()
    << " ms" << std::endl;
    
    return Y;
}
