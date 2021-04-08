//
//  Grid.hpp
//  CauchyProcess
//
//  Created by Rémi Carnec on 03/06/2020.
//  Copyright © 2020 Rémi Carnec. All rights reserved.
//

#ifndef Grid_hpp
#define Grid_hpp

#include "VectorHelper.hpp"

// This class represent a Grid of points. It is used to create a covering of a ball or a rectangle

enum TypeGrid
{
    Ball,
    Rectangle
};

class Grid {
    
private:
    TypeGrid typeGrid;
    double r;
    double delta;
    int dimension;
    int nbDomains;
    int nAxis;
    std::vector<double> sidesRect;
    std::vector<double> sidesDomains;
    std::vector<std::vector<double>> grid;
    std::vector<bool> activatedDomains;
    std::vector<std::vector<int>> neighboursIndices;
    
    
public:
    // Constructors for Ball and Rectangle
    Grid(double r, double delta, int dimension);
    Grid(std::vector<double> sidesRect, std::vector<double> sidesCover, int dimension);
    
    // Accessors
    TypeGrid getTypeGrid();
    double getRadius();
    double getDelta();
    int getDimension();
    int getNumDomains();
    std::vector<std::vector<double>> getGrid();
    std::vector<double> getSidesRect();
    std::vector<double> getSidesDomains();
    std::vector<bool> getActivatedDomains();
    std::vector<std::vector<int>> getNeighboursIndices();
    
    // Methods
    double computeVolume();
    int findIndex(std::vector<double> x);
    int findNbDomains(const std::vector<double>& x);
    int findNearestDomain(const std::vector<double>& x);
    std::vector<int> findIndices(const std::vector<double>& x);

};

#endif /* Grid_hpp */
