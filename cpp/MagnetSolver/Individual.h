#pragma once

#include <Eigen/Dense>
#include <vector>
#include "World.h"

// This is an individual, a component of the genetic algorithm. An individual is a single member of a population.

class Individual {

//private attributes
private:
	Eigen::Matrix<float, Eigen::Dynamic, 2> _Angles;											// A matrix with n rows and 2 columns to hold the angles of the tentacle.
	float _fitness;																				// The fitness of this particular solution.
																								
// No public Attriubutes.

// public methods
public:
	Individual(World world);																	// Constructor, uses the World object.
	~Individual();																				// Destructor

	// Accessors
	float& getFitness();																		// Accessor to return a reference to the fitness
	Eigen::Matrix<float, Eigen::Dynamic, 2>& getAngles();										// Accessor to retrun a reference to the angles

// private methods
private:
	void Randomise();																			// Randomises the joint angles.
	void DetermineFitness(World world);															// Measures the fitness of this soultion.
};