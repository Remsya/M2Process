//
//  PoissonPoint.cpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 18/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#include "PoissonPoint.hpp"

// Constructor
PoissonPoint::PoissonPoint(){}

PoissonPoint::PoissonPoint(double t, std::vector<double>& s, Grid *grid)
{
    T = t;
    S = s;
    
    // Initialize domain of influence
    if (grid->getTypeGrid() == TypeGrid::Rectangle)
    {
        // First corner of rectangle
        std::vector<double> vect = grid->getSidesRect();
        std::transform(vect.begin(), vect.end(), vect.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1./2));
        rectInfluence.push_back(vect);
        
        // Second corner
        std::transform(vect.begin(), vect.end(), vect.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
        rectInfluence.push_back(vect);
    }
    else if (grid->getTypeGrid() == TypeGrid::Ball)
    {
        std::vector<double> vect1(grid->getDimension(), -grid->getRadius());
        std::vector<double> vect2(grid->getDimension(), grid->getRadius());
        rectInfluence = {vect1, vect2};
    }
}

PoissonPoint::PoissonPoint(double t, std::vector<double>& s, std::vector<double>& corner1, std::vector<double>& corner2)
{
    // Check dimensions
    if (s.size() != corner1.size() || s.size() != corner2.size())
    {
        throw("Cannot create Poisson point: domain of influence is incoherent with dimension of S");
    }
    
    T = t;
    S = s;
    rectInfluence = {corner1, corner2};
}


// See if a point has influence relatively to a prior point
bool PoissonPoint::hasInfluence(const PoissonPoint& prevPoint, const StormFunction& shape)
{
    if (T < prevPoint.getT())
    {
        return true;
    }
    if (shape.getShape() == nameShapeFunction::Cauchy)
    {
        double lambda_1 = pow(T,1/shape.getAlpha());
        double lambda_0 = pow(prevPoint.getT(), 1/shape.getAlpha());
        double rho = lambda_1*lambda_0/pow(lambda_1 - lambda_0, 2)*pow(distance(S, prevPoint.getS()),2) - pow(shape.getScaleFactor(),2);
        
        // 1st Criteria for influence
        if (rho > 0)
        {
            // Center of ball of influence
            std::vector<double> center = linearCombination(S, prevPoint.getS(), lambda_1/(lambda_1 - lambda_0), -lambda_0/(lambda_1-lambda_0));
            
            return BallRectIntersect(center, sqrt(rho), shape);
        }
        else
        {
            return false;
        }
        
    }
    else if (shape.getShape() == nameShapeFunction::Laplace && shape.getAlpha() <= 2)
    {
        double alpha = shape.getAlpha();

        // Plane delimiting the real domains of influence P = {X \in R^d, X.n + D = 0}
        std::vector<double> n = substractVector(S, prevPoint.getS());
        double c = pow(shape.getScaleFactor()/norm(n),shape.getAlpha())*log(T/prevPoint.getT());
        double D;
        
        if (shape.getAlpha() == 2)
        {
            D = 1.0/2*(-pow(shape.getScaleFactor(),2)*log(T/prevPoint.getT()) + pow(norm(prevPoint.getS()),2) - pow(norm(S),2));
            
            return HalfSpaceRectIntersect(n, D, shape);
        }
        else if (shape.getAlpha() > 1)
        {
            double lambda;
            if (c == 1)
            {
                lambda = 1;
            }
            else if (c < 1)
            {
                lambda = pow(pow(2,-alpha)+c/2, 1/alpha);
            }
            else
            {
                lambda = (c-pow(2,alpha)*(1-alpha)-2*alpha+1)/(alpha*(pow(2,alpha-1)-1));
            }
            D = - lambda * pow(norm(n),2) - std::inner_product(n.begin(), n.end(), prevPoint.getS().begin(), 0.);
            
            return HalfSpaceRectIntersect(n, D, shape);
        }
        else if (alpha == 1)
        {
            double lambda;
            if (c >= 1)
            {
                return false;
            }
            else
            {
                lambda = (c+1)/2;
            }
            D = - lambda * pow(norm(n),2) - std::inner_product(n.begin(), n.end(), prevPoint.getS().begin(), 0.);
            
            return HalfSpaceRectIntersect(n, D, shape);
        }
        else
        {
            double lambda;
            if (c == 1)
            {
                lambda = 1;
            }
            else if (c < 1)
            {
                lambda = 1 - pow(pow(2,-alpha)-c/2,1/alpha); // Can add lambda_2 for sphere
            }
            else
            {
                return false;
            }
            D = - lambda * pow(norm(n),2) - std::inner_product(n.begin(), n.end(), prevPoint.getS().begin(), 0.);
            
            return HalfSpaceRectIntersect(n, D, shape);
        }
    }
    return true;
}

bool PoissonPoint::BallRectIntersect(std::vector<double>& center, double r, const StormFunction& shape)
{
    // Rectangle of influence
    std::vector<std::vector<double>> rect = {addConstant(center, -r), addConstant(center, r)};
    
    // 2nd Criteria for influence: rectangle must intersect
    for (int i = 0; i < shape.getDimension(); ++i)
    {
        if (rectInfluence[0][i] < rect[1][i] && rect[0][i] < rectInfluence[1][i])
        {
            rectInfluence[0][i] = std::max(rectInfluence[0][i], rect[0][i]);
            rectInfluence[1][i] = std::min(rectInfluence[1][i], rect[1][i]);
        }
        else
        {
            return false;
        }
    }
    
    return true;
}

bool PoissonPoint::HalfSpaceRectIntersect(std::vector<double>& n, double& D, const StormFunction& shape)
{
    // Check if all vertices S_i are on the same side of the plane
    std::vector<double> S_i = rectInfluence[0];
    double val_0 = std::inner_product(n.begin(), n.end(), S_i.begin(), 0.) + D > 0? 1 : -1, val = val_0;
    
    int i = 1;
    while (i < pow(2,shape.getDimension()) && val == val_0)
    {
        // Build i^th vertex
        for (int j = 0; j < shape.getDimension(); ++j)
        {
            // Index of j^th coordinates of i^th vertex, either 0 or 1
            int index = int(i/pow(2,j))%2;
            S_i[j] = rectInfluence[index][j];
        }
        
        // See if i^th vertex is on the same side as the previous vertices
        val = std::inner_product(n.begin(), n.end(), S_i.begin(), 0.) + D > 0? 1 : -1;
        ++i;
    }
    
    if (val == val_0 && val_0 == 1)
    {
        // Nothing to do, all rectangle is in the domain of influence
        return true;
    }
    else if (val == val_0 && val_0 == -1)
    {
        // The rectangle has no intersection with the half space of influence, thus the new Poisson point has no influence anymore
        return false;
    }
    else
    {
        // The rectangle intersect with the plane, build new rectangle coordinate by coordinate
        std::vector<std::vector<double>> newRectInfluence = rectInfluence;
        for (int k = 0; k < shape.getDimension(); ++k)
        {
            // Initialize
            if (n[k] >= 0)
            {
                newRectInfluence[0][k] = rectInfluence[1][k];
            }
            else
            {
                newRectInfluence[1][k] = rectInfluence[0][k];
            }
            
            // Update coordinates looking at all edges
            for (int i = 0; i < pow(2, shape.getDimension() - 1); ++i)
            {
                // Build edge i
                std::vector<double> edge_i(shape.getDimension(), 0.);
                for (int j = 0; j < shape.getDimension(); ++j)
                {
                    if (j < k)
                    {
                        int index = int(i/pow(2,j))%2;
                        edge_i[j] = rectInfluence[index][j];
                    }
                    else if (j > k)
                    {
                        int index = int(i/pow(2,j-1))%2;
                        edge_i[j] = rectInfluence[index][j];
                    }
                    
                }
                
                // Compute intersection between edge i and plane
                double intersect = intersectCoordinate(n, D, edge_i, k);
                
                // Update newRectangle
                newRectInfluence[0][k] = std::min(newRectInfluence[0][k], std::max(rectInfluence[0][k], intersect));
                newRectInfluence[1][k] = std::max(newRectInfluence[1][k], std::min(rectInfluence[1][k], intersect));
            }
            rectInfluence = newRectInfluence;
        }
    }
    
    return true;
}

// Check if x is in domain of influence of a Poisson point
bool PoissonPoint::contributes(const std::vector<double>& x) const
{
    bool cont = true;
    for (int i = 0; cont && i < x.size() ; ++i)
    {
        cont = (x[i] >= rectInfluence[0][i] && x[i] <= rectInfluence[1][i]);
    }
    return cont;
}

// Deal with poisson points for data processing
void extractSValues(std::vector<PoissonPoint>& x, std::vector<std::vector<double>>& S)
{
    // Explore poisson points
    for (PoissonPoint& p: x)
    {
        S.push_back(p.getS());
    }
}

void extractTValues(std::vector<PoissonPoint>& x, std::vector<double>& T)
{
    // Explore poisson points
    for (PoissonPoint& p: x)
    {
        T.push_back(p.getT());
    }
}

void extractDomainsInfluence(std::vector<PoissonPoint>& x, std::vector<std::vector<double>>& rect)
{
    // Explore poisson points
    for (PoissonPoint& p: x)
    {
        rect.push_back(p.getRectInfluence()[0]);
        rect.push_back(p.getRectInfluence()[1]);
    }
}


// Find intersection between plane P = {X \in R^d, X.n + D = 0} and the line led by edge "edgeVertex"
double intersectCoordinate(std::vector<double>& n, double& D, std::vector<double>& edgeVertex, int index)
{
    // Check
    if (index < 0 || index >= edgeVertex.size())
    {
        throw("Cannot compute intersection of plane and line. Index");
    }
    
    // Compute intersection
    return (- D - std::inner_product(n.begin(), n.end(), edgeVertex.begin(), 0.) + n[index] * edgeVertex[index])/n[index];
}
