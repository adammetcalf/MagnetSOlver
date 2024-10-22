// Individual class. An individual is a single member of a population


#include "Individual.h"


// Constructor
Individual::Individual(World world) : world_(world) {

	Individual::Randomise();
	Individual::DetermineFitness();
}


// Destructor
Individual::~Individual() {

}

// Update the field and recalculate tentacle/ soultion fitness
void Individual::UpdateField(double Bx, double  Bz) {

	Bx_ = Bx;
	Bz_ = Bz;
	Individual::DetermineFitness();
}

// Accessor to return (a reference to) the fitness
double& Individual::GetFitness(){

	return fitness_;
}

void Individual::Randomise() {

	// Seed random number generator
	std::srand(std::time(0));

	// define the limits
	double lowerLimit = -0.0018;
	double upperLimit = 0.0018;

	// Randomize a and b between the limits
	Bx_ = lowerLimit + (upperLimit - lowerLimit) * (std::rand() / (double)RAND_MAX);
	Bz_ = lowerLimit + (upperLimit - lowerLimit) * (std::rand() / (double)RAND_MAX);
}

// Method to determine the fitness fo this solution
void Individual::DetermineFitness() {

	// Update the field
	world_.UpdateField(Bx_, Bz_);

	// get the desired and actual positions
	const auto& CurrentX = world_.GetCurrentX();
	const auto& CurrentZ = world_.GetCurrentZ();
	const auto& DesX = world_.GetDesiredX();
	const auto& DesZ = world_.GetDesiredZ();

	// Calculate difference between current and desired positions
	Eigen::VectorXd diffX = CurrentX - DesX;
	Eigen::VectorXd diffZ = CurrentZ - DesZ;

	// Concatenate X and Z differences
	Eigen::VectorXd diff;
	diff.resize(diffX.size() + diffZ.size());
	diff << diffX, diffZ;

	// Compute the norm of the concatenated difference vector and square it --- vector norm
	fitness_ =  std::pow(diff.norm(), 2);
}

// Accessor to return the fields
std::vector<double> Individual::GetField() {

	std::vector<double> Fields;

	Fields[0] = Bx_;
	Fields[1] = Bz_;

	return Fields;
}
