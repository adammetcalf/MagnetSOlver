// This is a joint, which is represented by a single standard Denavit Hartenberg frame

#define _USE_MATH_DEFINES																									// To ensure that M_PI works
#include "Joint.h"
#include <cmath>																											// For M_PI


// Constructor
Joint::Joint(float theta, float alpha, float a, float d)
	: _theta(deg2rad(theta)), _alpha(deg2rad(alpha)), _a(a), _d(d), _Frame(Eigen::Matrix4f::Zero()) {						// Member initialisation list - inits the class members. Note angles are immediately converted to rad
	
	Joint::evaluateFrame();
}

// Defualt constructor																									
Joint::Joint()																												// COnstructor overloading - this object is instantiated when no iput arguements are provided
	: _theta(0.0f), _alpha(0.0f), _a(0.0f), _d(0.0f), _Frame(Eigen::Matrix4f::Identity()) {
	evaluateFrame();
}

// Destructor
Joint::~Joint() {
	// no action necessary.
};


// Return a reference to the DH Frame
Eigen::Matrix4f& Joint::getFrame() {

	return _Frame;
};


// Evaluate the DH Frame
void Joint::evaluateFrame() {

	// overwrite frame with newly calculated parameters.
	_Frame << cos(_theta), -sin(_theta) * cos(_alpha), sin(_theta)* sin(_alpha), _a* cos(_theta),
		sin(_theta), cos(_theta)* cos(_alpha), -cos(_theta) * sin(_alpha), _a* sin(_theta),
		0, sin(_alpha), cos(_alpha), _d,
		0, 0, 0, 1;

	// Set any tiny elements to zero
	_Frame = _Frame.unaryExpr([](float elem) { return fabs(elem) < 1.0e-7 ? 0.0f : elem; });

	// Eigen::Matrix.unaryExpr applies an operation to each element in the array.
	// The operation is a lambda
	// The lambda checks whether the absoulte value is less than small number, and returns 0.0 if so. if not, it returns the element
}


// Update the angles
void Joint::UpdateAngles(std::vector<float> Angles) {

	_theta = deg2rad(Angles[0]);
	_alpha = deg2rad(Angles[1]);

	Joint::evaluateFrame();																									//revaluate frame after joint angle update.
}


// Convert degrees to radians
float Joint::deg2rad(float degrees) {
	return degrees * M_PI / 180.0f;
}