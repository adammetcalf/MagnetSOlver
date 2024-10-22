// This class contains a parent selected from the population for breeding

#pragma once

#include <vector>
#include <cstdlib>										// for the randomise
#include <ctime>										// for the randomise
#include "Individual.h"


class Parent {

// Private propeties
private:
	Individual progenator_;

// Public Methods
public:

	// Constructor
	Parent(std::vector<Individual>& popArray); // takes in a reference to the vector of individuals (population array)

	// Destructor
	~Parent();

	// return the field of the progenator, which is actually the only thing we care about solution-wise
	std::vector<double> GetField();

// Private Methods
private:

	// Obtain a parent from the population array
	void GetParent(std::vector<Individual>& popArray);

	// biased random function - select an individual using biased random mehtod
	void BiasedRandom(std::vector<Individual>& popArray);

	// tournament function - randomly selects 2 individuals and chooses the fittest
	void Tournament(std::vector<Individual>& popArray);

};