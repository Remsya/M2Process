//
//  PoissonPoint.hpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 18/05/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#ifndef PoissonPoint_hpp
#define PoissonPoint_hpp

#include "StormFunction.hpp"
#include "Grid.hpp"
#include "M2ProcessHelper.hpp"

// This class represent a Poisson Process Point (T,S), it is used to eliminate useless candidates for continuous Max Stable simulation

class PoissonPoint {
    
private:
    double T;
    std::vector<double> S;
    std::vector<std::vector<double>> rectInfluence;
    
public:
    // Default constructor
    PoissonPoint();
    
    // Constructors
    PoissonPoint(double t, std::vector<double>& s, Grid *grid);
    PoissonPoint(double t, std::vector<double>& s, std::vector<double>& corner1, std::vector<double>& corner2);
    
    // Accessors
    inline const double& getT() const {return T;}
    inline const std::vector<double>& getS() const {return S;}
    inline const std::vector<std::vector<double>>& getRectInfluence() const {return rectInfluence;}
    
    // Method: See where the Poisson point may have influence
    bool hasInfluence(const PoissonPoint& prevPoint, const StormFunction& shape);
    bool BallRectIntersect(std::vector<double>& center, double r, const StormFunction& shape);
    bool HalfSpaceRectIntersect(std::vector<double>& n, double& D, const StormFunction& shape);
    bool contributes(const std::vector<double>& x) const;
};

// Intersection of domains of influence
double intersectCoordinate(std::vector<double>& n, double& D, std::vector<double>& edgeVertex, int index);

// Extract values
void extractSValues(std::vector<PoissonPoint>& x, std::vector<std::vector<double>>& S);
void extractTValues(std::vector<PoissonPoint>& x, std::vector<double>& T);
void extractDomainsInfluence(std::vector<PoissonPoint>& x, std::vector<std::vector<double>>& rect);

#endif /* PoissonPoint_hpp */
