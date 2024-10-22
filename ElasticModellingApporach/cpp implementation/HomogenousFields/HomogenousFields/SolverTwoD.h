// This class numerically solves for the position of a flexible
// magnetic tentacle under homogenous fields


#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <Eigen/Dense>


class SolverTwoD {

// Private Properties 
private:

    int n_;                                // Number of nodes in the tentacle
    double L_;                             // Total length of the tentacle (m)
    double Bxin_;                          // X-component of the magnetic field (T)
    double Bzin_;                          // Z-component of the magnetic field (T)
    Eigen::MatrixXd magMoments_;           // Magnetic moments at each node (n x 2 matrix)
    int orientation_;                      // Starting orientation (0 for horizontal, 1 for vertical downwards)
    Eigen::VectorXd xpos_;                 // X positions of the tentacle nodes
    Eigen::VectorXd zpos_;                 // Z positions of the tentacle nodes 
    Eigen::VectorXd theta_;                // Rotation angle at each node
    Eigen::VectorXd theta_p_;              // First derivative of theta
    Eigen::VectorXd theta_p_old_;          // Previous iteration of theta_p_
    Eigen::MatrixXd r_;                    // Position vectors at each node (n x 2 matrix)
    Eigen::MatrixXd m_def_;                // Deformed magnetic moment vectors
    Eigen::MatrixXd m_hat_;                // Derivative of m_def_ with respect to theta
    double const1_;                        // Precomputed constant: (rho * g * A) / (E * I)
    double const2_;                        // Precomputed constant: 1 / (E * I)
    double const3_;                        // Precomputed constant: (E * I)
    double delta_s_;                       // Link length between nodes
    bool use_external_forces_;             // Flag to indicate if external forces are used
    const Eigen::VectorXd Fx_;             // External forces in the x-direction at each node
    const Eigen::VectorXd Fz_;             // External forces in the z-direction at each node

    // Private methods
    void initialise();                     // Initializes variables and precomputes constants
    void solve();                          // Iteratively solves for the tentacle positions

    //cumulutive sum helper function
    Eigen::VectorXd cumulativeSum(const Eigen::VectorXd&);

public:

    // Constructor for intial solve
    SolverTwoD(int n, double L, double Bxin, double Bzin, const Eigen::MatrixXd& magMoments, int orientation);

    // overloaded operater to accept external force inputs also. This is used for the admittance control
    SolverTwoD(int n, double L, double Bxin, double Bzin, const Eigen::MatrixXd& magMoments, int orientation,
        const Eigen::VectorXd& Fx, const Eigen::VectorXd& Fz);
    
    // Destructor
    ~SolverTwoD();

    // return the xpositions of the nodes
    const Eigen::VectorXd& getXPositions() const;

    // return the zpositions of the nodes
    const Eigen::VectorXd& getZPositions() const;
};