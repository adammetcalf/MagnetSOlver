// The world contains the magnetically actuated elastica and the actuating magnetic field.
// The world has a field and a tentacle solver (with all the necessary supporting properties).
// The world throws an updated field at the solver, and returns the positions

#include "World.h"

// Constructor - No external forces
World::World(int n, double L, const Eigen::MatrixXd& magMoments, int orientation, const Eigen::VectorXd& DesX, const Eigen::VectorXd& DesZ) :
	n_(n), L_(L), magMoments_(magMoments), orientation_(orientation), DesX_(DesX), DesZ_(DesZ), Bx_(0.0), Bz_(0.0), Fx_(Eigen::VectorXd::Zero(1)), Fz_(Eigen::VectorXd::Zero(1)), ExternalForces_(false) {

}

// Constructor overloaded with external forces
World::World(int n, double L, const Eigen::MatrixXd& magMoments, int orientation, const Eigen::VectorXd& DesX, const Eigen::VectorXd& DesZ, const Eigen::VectorXd& Fx, const Eigen::VectorXd& Fz) :
	n_(n), L_(L), magMoments_(magMoments), orientation_(orientation), DesX_(DesX), DesZ_(DesZ), Bx_(0.0), Bz_(0.0), Fx_(Fx), Fz_(Fz), ExternalForces_(true){

}


// Destructor
World::~World() {

}

// Update Field method -- update Bx and Bz, then re-evaluate the tentacle.
void World::UpdateField(double Bx, double Bz) {

	Bx_ = Bx;
	Bz_ = Bz;

	World::EvaluateWorld();
}

// Evaluate the model using the conditions specific to this world
void World::EvaluateWorld() {

	// If using external forces (Admittance Control)
	if (ExternalForces_){

		// instantiate solverTwoD class
		SolverTwoD model(n_, L_, Bx_, Bz_, magMoments_, orientation_, Fx_, Fz_);
		CurrentX_ = model.getXPositions();
		CurrentZ_ = model.getZPositions();
	}
	// No external forces present
	else {

		// instantiate solverTwoD class
		SolverTwoD model(n_, L_, Bx_, Bz_, magMoments_, orientation_);
		CurrentX_ = model.getXPositions();
		CurrentZ_ = model.getZPositions();
	}
}

// Accessor to get a reference to CurrentX
const Eigen::VectorXd& World::GetCurrentX() {

	return CurrentX_;
}

// Accessor to get a reference to CurrentY
const Eigen::VectorXd& World::GetCurrentZ() {

	return CurrentZ_;
}

// Accessor to get a reference to DesiredX
const Eigen::VectorXd& World::GetDesiredX() {

	return DesX_;
}

// Accessor to get a reference to DesiredZ
const Eigen::VectorXd& World::GetDesiredZ() {

	return DesZ_;
}
