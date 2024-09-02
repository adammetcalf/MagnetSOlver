#pragma once

#include <Eigen/Dense>
#include <vector>


// This is a joint, which is represented by a single standard Denavit Hartenberg frame

class Joint {

// Private Properties
private:
	float _a;																		// DH Parameter (in m)
	float _d;																		// DH Parameter (in m)
	float _theta;																	// DH Parameter (in radians)
	float _alpha;																	// DH Parameter (in radians)
	Eigen::Matrix4f _Frame;															// DH Frame

// Public Properties: None.

// Public Methods
public:

	Joint(float theta, float alpha, float a, float d);								// Constructor
	Joint();																		// Default constrcutor
	~Joint();																		// Destructor
	void UpdateAngles(std::vector<float> Angles);									// Update the internally held angles of the joint
	Eigen::Matrix4f& getFrame();													// Accessor to return a reference to the frame (No copies!!!!)									

// Private methods
private:

	void evaluateFrame();															// Private function to evaluate the frame
	float deg2rad(float degrees);													// Helper function to convert degrees to radians
};



