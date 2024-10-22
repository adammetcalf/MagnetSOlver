// Individual class. An individual is a single member of a population

#pragma once
#include <Eigen/Dense>
#include <vector>
#include <tuple>
#include <cmath>										// For pow()
#include <cstdlib>										// for the randomise
#include <ctime>										// for the randomise
#include "World.h"


class Individual {

	// private properties
private:
	double fitness_ = 10000.0;							// a fitness measure of this particular solution
	double Bx_ = 0.0;									// x component of the magnetic field
	double Bz_ = 0.0;									// y component of the magnetic field
	World world_;										// An instance of the world class contained by this object.

// public methods
public:

	// Constructor
	Individual(World world);

	// Destructor
	~Individual();

	// Method to update the magnetic field
	void UpdateField(double Bx, double Bz);

	// Accessor to return (a reference to) the fitness
	double& GetFitness();

	// Accessor to return the fields
	std::vector<double> GetField();

	//private methods
private:

	// Method to randomise the field
	void Randomise();

	// Method to determine the fitness fo this solution
	void DetermineFitness();
};
