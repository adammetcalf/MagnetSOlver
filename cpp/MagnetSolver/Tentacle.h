#pragma once

#include <vector>
#include "Joint.h"																				// To enable a vector of Joint objects.
#include <Eigen/Dense>

// This is a tentacle class, composed of Joints and magnetic moments.


class Tentacle {

// private attributes
private:

	std::vector<Joint> _Joints;																	// A vector to hold the joint objects.
	Eigen::Matrix<float, Eigen::Dynamic, 2> _Angles;											// A matrix with n rows and 2 columns to hold the angles.
	float _LinkLength;																			// the distance between each joint.
	std::vector<Eigen::Matrix4f> _HGM;															// A vector of 4x4 matrices (Homogeneous transformation matrices)
	//Eigen::Matrix<double, 4, 4, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic> _HGM;		// A 4x4 matrix x n depth, essentially an array of homogeneous transformation matrices.

// No public properties


//public methods
public:

	Tentacle(float LinkLength, Eigen::Matrix<float, Eigen::Dynamic, 2> Angles);					// Constructor
	~Tentacle();																				// Destructor
	void UpdateAngles(Eigen::Matrix<float, Eigen::Dynamic, 2> newAngles);						// Function to update the joint angles of the tentacle.
	std::vector<Eigen::Matrix4f>& GetHGMs();													// Accessor to get a reference to the homogeneous transformation matrices.
	Eigen::Matrix<float, Eigen::Dynamic, 2>& GetAngles();										// Accessor to get a reference to the joint angles.

// private methods
private:

	void EvaluateHGMs();																		// Method to evaluate the homogeneous transformation matrices

};