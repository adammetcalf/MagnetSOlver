// The world contains the magnetically actuated elastica and the actuating magnetic field.
// The world has a field and a tentacle solver (with all the necessary supporting properties).
// The world throws an updated field at the solver, and returns the positions

#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>
#include "SolverTwoD.h"

class World {

// private properties
private:
    
    Eigen::VectorXd CurrentX_;             // Current X positions of the tentacle (m)
    Eigen::VectorXd CurrentZ_;             // Current Z positions of the tentacle (m)
    const Eigen::VectorXd DesX_;           // Desired X positions of the tentacle (m)
    const Eigen::VectorXd DesZ_;           // Desired Z positions of the tentacle (m)
    
    // tentacle properties
    int orientation_;                      // Starting orientation (0 for horizontal, 1 for vertical downwards)
    double L_;                             // Total length of the tentacle (m)
    const Eigen::MatrixXd magMoments_;     // Magnetic moments at each node (n x 2 matrix)
    int n_;                                // Number of nodes in the tentacle

    // Mag proeprties
    double Bx_;                             // Field in x direction (T)
    double Bz_;                             // Field in z-direction (T)

    // External forces
    Eigen::VectorXd Fx_;                   // External forces in the x-direction at each node
    Eigen::VectorXd Fz_;                   // External forces in the z-direction at each node
    bool ExternalForces_;                  // Flag to determine whether to use external forces in the mdoel.


// public methods
public:

    // constructor
    World(int n, double L, const Eigen::MatrixXd& magMoments, int orientation, const Eigen::VectorXd& DesX, const Eigen::VectorXd& DesZ);

    // constructor overloaded with external forces
    World(int n, double L, const Eigen::MatrixXd& magMoments, int orientation, const Eigen::VectorXd& DesX, const Eigen::VectorXd& DesZ, const Eigen::VectorXd& Fx, const Eigen::VectorXd& Fz);

    //Destructor
    ~World();

    // Update Field
    void UpdateField(double Bx, double Bz);


    // Accessor to get a reference to CurrentX
    const Eigen::VectorXd& GetCurrentX();

    // Accessor to get a reference to CurrentY
    const Eigen::VectorXd& GetCurrentZ();

    // Accessor to get a reference to DesiredX
    const Eigen::VectorXd& GetDesiredX();

    // Accessor to get a reference to DesiredZ
    const Eigen::VectorXd& GetDesiredZ();

// private
private:

    // Evaluate World
    void EvaluateWorld();
};