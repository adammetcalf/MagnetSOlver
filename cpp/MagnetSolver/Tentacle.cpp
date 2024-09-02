// This is the tentacle, which encompasses n joints.

#define _USE_MATH_DEFINES																		// To ensure that M_PI works
#include "Tentacle.h"
#include <cmath>


// Constructor
Tentacle::Tentacle(float LinkLength, Eigen::Matrix<float, Eigen::Dynamic, 2> Angles) 
	: _LinkLength(LinkLength), _Angles(Angles) {												// Member initialisation list - inits the class members.

    _Joints.resize(Angles.rows());                                                              // Resize the Joints vector to hold the same number of joints as rows in the Angles matrix

    for (int i = 0; i < Angles.rows(); ++i) {                                                   // Create 1 joint for each row in the Angles Matrix
        if (i == 0) {
            
            _Joints[i] = Joint(Angles(i, 0), Angles(i, 1), 0.0f, 0.0f);                          // The 1st joint at (0,0,0) has no length displacement
        }
        else {
            
            _Joints[i] = Joint(Angles(i, 0), Angles(i, 1), 0.0f, LinkLength);                    // Every other joint is displaced LinkLenght along DH parameter 'd'
        }
    }

    Tentacle::EvaluateHGMs();                                                                   // Evaluate the Homogeneous Transformation Matrices      
};

// Destructor
Tentacle::~Tentacle() {                                                                         // Nothing to do on clearup.
	
};

// Accessor to get the a refrenece to the HGMs
std::vector<Eigen::Matrix4f>& Tentacle::GetHGMs() {

	return _HGM;
};

// Accessor to get a reference to the joint angles
Eigen::Matrix<float, Eigen::Dynamic, 2>& Tentacle::GetAngles() {

	return _Angles;
}

// Evaluate Homogeneous Transformation Matrices
void Tentacle::EvaluateHGMs() {
    Eigen::Matrix4f Origin = Eigen::Matrix4f::Identity();                                       // Define using a 4x4 identity matrix
    _HGM.clear();                                                                               // Clear the vector (to ensure it's empty before we start)
    _HGM.reserve(_Joints.size());                                                               // Reserve vector length memory for efficiency

    Eigen::Matrix4f HGM = Origin;                                                               // Initialize the first HGM

    for (int i = 0; i < _Joints.size(); ++i) {
        
        Joint& J = _Joints[i];                                                                  // Get a reference to the ith joint
        Eigen::Matrix4f Frame = J.getFrame();                                                   // Get the joint frame

        if (i == 0) {
            HGM = Origin * Frame;                                                               // Initialize using transform from origin to frame 1
        }
        else {
            HGM = HGM * Frame;                                                                  // Transform from previous HGM to next HGM
        }
        _HGM.push_back(HGM);                                                                    // Store the result in the vector of HGMs
    }
}


// Update the joint angles
void Tentacle::UpdateAngles(Eigen::Matrix<float, Eigen::Dynamic, 2> newAngles) {

    this->_Angles = newAngles;
}