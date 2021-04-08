//
//  M2ProcessHelper.cpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 21/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include "M2ProcessHelper.hpp"


// Beta function: Beta(x,y) = Gamma(x)*Gamma(y)/Gamma(x+y)
double beta(double x, double y)
{
    return tgamma(x)*tgamma(y)/tgamma(x+y);
}

// Returns the binomial coefficient C(k,n) using Pascal triangle
double comb(int k, int n)
{
    if (k == 0 || k==n)
    {
        return 1;
    }
    else if (k > n/2.)
    {
        return comb(n-k,n);
    }
    else
    {
        return comb(k-1,n-1)+comb(k,n-1);
    }
}

// Compute the volume D1 = {(s,t) \in B_r\R X R_+^*, t < C_N*f(d(s,B_r))}, for rectangle grid
double computeVolumeD1(const StormFunction& shape, Grid *grid, double C_N)
{
    // Function only for rectangle type grid
    if (grid->getTypeGrid() != TypeGrid::Rectangle)
    {
        throw("Function computeVolumeD1 only for Rectangle window");
    }
    
    // Volume of unit d-ball
    int d = shape.getDimension();
    double w_d = pow(M_PI, d/2.)/tgamma(d/2.+1);
    
    // Compute volume of D1
    std::vector<double> sidesRect = grid->getSidesRect();
    double prodLengths = std::accumulate(sidesRect.begin(), sidesRect.end(), 1., std::multiplies<double>());
    double vol = C_N*(w_d*pow(grid->getRadius(),d) - prodLengths);
    
    return vol;
}

// Compute the volume D2 = {(s,t) \in B_r^c X R_+^*, t < C_N*f(d(s,B_r))}
double computeVolumeD2(const StormFunction& shape, Grid *grid, double C_N)
{
    // Compute volume
    int d = shape.getDimension();
    double w_d = pow(M_PI, d/2.)/tgamma(d/2.+1);
    double r = grid->getRadius();
    double vol = 0;
    for (int i = 0; i < d; ++i)
    {
        vol += C_N*d*w_d*comb(i,d-1)*pow(r,d-1-i)*shape.getMoment(i);
    }
    
    return vol;
}

// Compute p_I*Volume
std::vector<double> weightsAndVolume(const StormFunction& shape, Grid* grid, double& volume, double C)
{
    std::vector<double> L = grid->getSidesRect();
    std::vector<double> p_I(pow(2, grid->getDimension()), 0.);
    std::vector<double> w(grid->getDimension()+1, 0);
    
    // Build unit volumes
    for (int i = 1; i < grid->getDimension()+1; ++i)
    {
        w[i] = pow(M_PI, i/2.)/tgamma(i/2.+1);
    }
    
    // Build p_i
    for (int i = 1; i < p_I.size(); ++i)
    {
        int cardinal = 0;
        double prodLengths = 1;
        for (int j = 0; j < grid->getDimension(); ++j)
        {
            if (int(i/pow(2,j))%2 == 0)
            {
                prodLengths *= L[j];
            }
            else
            {
                ++cardinal;
            }
        }
        p_I[i] = prodLengths*cardinal*shape.getMoment(cardinal-1)*w[cardinal];
        volume += p_I[i];
    }
    volume *= C;
    
    return p_I;
}

// Counts the number of "1" in the binary expression of n
int countOnes (int n)
{
    int count=0;
    while (n!=0)
    {
        n = n & (n-1);
        count++;
    }
    return count;
}

double distanceFromGrid(std::vector<double>& x, Grid *grid)
{
    if (x.size() != grid->getDimension())
    {
        throw("x size must be equal to dimension");
    }
    
    std::vector<double> rectangle = grid->getSidesRect();
    double distance = 0;
    
    // Compute square of distance
    for (int i = 0; i < x.size(); ++i)
    {
        distance += abs(x[i]) > rectangle[i]/2? pow(abs(x[i])-rectangle[i]/2,2) : 0;
    }
    
    return sqrt(distance);
}



// Simulation functions

// Simulate a random variable that has a density g
std::vector<double> simulateDensityBall(std::default_random_engine& gen, const StormFunction& shape, double r)
{
    // Parameters
    double a = shape.getScaleFactor();
    double alpha = shape.getAlpha();
    int n = shape.getDimension();
    
    // Need to compute the sum of the StormFunction moments
    double sum = shape.getCumMoment();
    
    // Need binomial and uniform distributions for the procedure to simulate k (to decide which f_k we use)
    std::binomial_distribution<>binomialDistribution(n-1,1/(r+1));
    std::uniform_real_distribution<>unifDistribution(0,1);
    
    // Initialize random Variables, simulate until m_X > U * (m_0 + ... + m_{d-1})
    int X = binomialDistribution(gen);
    double U = unifDistribution(gen);
    while (shape.getMoment(X)/sum < U)
    {
        X = binomialDistribution(gen);
        U = unifDistribution(gen);
    }

    // Simulate Y ~ r + f_X
    double Y;
    if (shape.getShape() == nameShapeFunction::Laplace)
    {
        // Generalized Gamma distribution
        Y = r + generalizedGammaDist(gen, a, X+1, alpha);
    }
    else if (shape.getShape() == nameShapeFunction::Cauchy)
    {
        // Simulate w
        std::gamma_distribution<>gammaDist(alpha-(X+1)/2.,1);
        double w = gammaDist(gen);
        
        // Generalized Gamma distribution with a' = a/sqrt(w) and alpha = 2
        Y = r + generalizedGammaDist(gen, a/sqrt(w), X+1, 2);
    }
    else
    {
        Y = r;
    }
    
    // Build coordinates of S {x_1,...,x_n} with radius Y
    
    // Build Multivariate normal distribution
    std::normal_distribution<double>normalDistribution(0,1);
    auto normalGenerator = std::bind(normalDistribution, std::ref(gen));
    
    // Build S
    std::vector<double> S(n);
    std::generate(S.begin(), S.end(), normalGenerator);
    
    // Multiply by Y/norm(S)
    double normalization = norm(S);
    std::transform(S.begin(), S.end(), S.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, Y/normalization));
    
    return S;
}

std::vector<double> simulateDensityRectangle(std::default_random_engine& gen, const StormFunction& shape, Grid *grid, const std::vector<double>& p_I)
{
    // Parameters
    int dim = shape.getDimension();
    std::vector<double> rectangle = grid->getSidesRect();
    
    // Simulate Subset I for k-faces
    std::discrete_distribution<> subsetIndexDistribution(p_I.begin(), p_I.end());
    int indexI = subsetIndexDistribution(gen);
    int cardI = countOnes(indexI);
    
    // Vector I = [1_{i \in I}, i \in {1,...,d} ]
    std::vector<int> I(dim);
    for (int j = 0; j < dim; ++j)
    {
        I[j] = int(indexI/pow(2,j))%2;
    }
    
    // Simulate distance from rectangle
    double distance = simulateFk(gen, shape, cardI - 1);
    
    // Simulate Y on cardI-sphere of radius 'distance'
    std::vector<double> Y = unifSphere(gen, cardI, distance);
    int counter = 0;
    
    // Build S
    std::uniform_real_distribution<> unifDistribution(0,1);
    std::vector<double> S(dim);
    for (int j = 0; j < dim; ++j)
    {
        if (I[j] == 1)
        {
            S[j] = Y[counter] + rectangle[j]/2 * (Y[counter] >= 0? 1 : -1);
            ++counter;
        }
        else
        {
            S[j] = rectangle[j]*(unifDistribution(gen)-1./2);
        }
    }
    
    return S;
}

double simulateFk(std::default_random_engine& gen, const StormFunction& shape, int k)
{
    if (shape.getShape() == nameShapeFunction::Laplace)
    {
        return generalizedGammaDist(gen, shape.getScaleFactor(), k+1, shape.getAlpha());
    }
    else if (shape.getShape() == nameShapeFunction::Cauchy)
    {
        std::gamma_distribution<>gammaDist(shape.getAlpha()-(k+1)/2.,1);
        double w = gammaDist(gen);
        return generalizedGammaDist(gen, shape.getScaleFactor()/sqrt(w), k+1, 2);
    }
    else
    {
        return 0;
    }
}

// Uniform distribution on a n-Sphere:
// If Y is uncorrelated multivariate normal distribution, then Y/|Y|*U(0,1)^{1/n} has uniform distribution on unit n-Ball
std::vector<double> unifSphere(std::default_random_engine& gen, int dimension, double r)
{
    // Need Uniform and normal random variable generators
    std::normal_distribution<double>normalDistribution(0,1);
    std::uniform_real_distribution<double>unifDistribution(0,1);

    // Normal random variable Generator
    auto normalGenerator = std::bind(normalDistribution, std::ref(gen));
    
    // Build S
    std::vector<double> S(dimension);
    std::generate(S.begin(), S.end(), normalGenerator);
    
    // Build uniform distribution in n-Ball with radius r
    std::transform(S.begin(), S.end(), S.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, r));
    
    return S;
}

// Simulate uniformDistribution in Ball or Rectangle
std::vector<double> unifGrid(std::default_random_engine& gen, Grid *grid)
{
    // Initialize S
    std::vector<double> S(grid->getDimension());
    
    if (grid->getTypeGrid() == TypeGrid::Ball)
    {
        // Need Uniform and normal random variable generators
        std::normal_distribution<double>normalDistribution(0,1);
        std::uniform_real_distribution<double>unifDistribution(0,1);

        // Normal random variable Generator
        auto normalGenerator = std::bind(normalDistribution, std::ref(gen));
        
        // Build S
        std::generate(S.begin(), S.end(), normalGenerator);
        
        // Build uniform distribution in n-Ball with radius r
        double radius = grid->getRadius() * pow(unifDistribution(gen),1./grid->getDimension()) / norm(S);
        std::transform(S.begin(), S.end(), S.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, radius));
    }
    else if (grid->getTypeGrid() == TypeGrid::Rectangle)
    {
        // Need Uniform and normal random variable generators
        std::uniform_real_distribution<double>unifDistribution(0,1);
        
        // Lengths (L_i) of rectangle
        std::vector<double> sidesRect = grid->getSidesRect();
        
        // Fill with U(-L_i/2,L_i/2)
        for (int i = 0; i < grid->getDimension(); ++i)
        {
            S[i] = sidesRect[i]*(unifDistribution(gen)-1./2);
        }
    }
    
    return S;
}

// If X~Generalized Gamma(a,beta,alpha) and W~Gamma(beta/alpha, a^alpha), then X~W^{1/alpha}
double generalizedGammaDist(std::default_random_engine& gen, double a, double beta, double alpha)
{
    // Simulate W~Gamma(beta/alpha, a^alpha)
    std::gamma_distribution<>gammaDist(beta/alpha, pow(a,alpha));
    double w = gammaDist(gen);
    
    // Return w^{1/alpha} to get the desired distribution
    return pow(w,1/alpha);
}
