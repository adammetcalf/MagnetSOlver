// This class contains a parent selected from the population for breeding

#include "Parent.h"


// Constructor
Parent::Parent(std::vector<Individual>& popArray) : progenator_(popArray[0]) {

	// Seed random number generator
	std::srand(std::time(0));

	// replace initilaised porgenator with actual useful entry
	Parent::GetParent(popArray);
}

// Destructor
Parent::~Parent() {

}

// return the field of the progenator, which is actually the only thing we care about solution-wise
std::vector<double> Parent::GetField() {

	return progenator_.GetField();
}

// Obtain a parent from the population array
void Parent::GetParent(std::vector<Individual>& popArray) {

	// get random float between 0 and 1
	float x = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);

	if (x > 0.5) {
		// use toournament
		Parent::Tournament(popArray);

	}
	else {
		// use biased random
		Parent::BiasedRandom(popArray);
	}
}

// biased random function - select an individual using biased random mehtod
void Parent::BiasedRandom(std::vector<Individual>& popArray) {

	// get random float between 0 and 1
	float bias = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);

	if (bias >= 0.7) { //30% chance
		progenator_ = popArray[0];
	}
	else if (bias >= 0.5) { //20% chance
		progenator_ = popArray[1];
	}
	else if (bias >= 0.35) { //15% chance
		progenator_ = popArray[2];
	}
	else if (bias >= 0.25) { //10% chance
		progenator_ = popArray[3];
	}
	else if (bias >= 0.2) { //5% chance
		progenator_ = popArray[4];
	}
	else if (bias >= 0.15) { //5% chance
		progenator_ = popArray[5];
	}
	else if (bias >= 0.1) { //5% chance
		progenator_ = popArray[6];
	}
	else if (bias >= 0.05) { //5% chance
		progenator_ = popArray[7];
	}
	else if (bias >= 0.025) { //2.5% chance
		progenator_ = popArray[8];
	}
	else { //2.5% chance
		progenator_ = popArray[9];
	}

}

// tournament function - randomly selects 2 individuals and chooses the fittest
void Parent::Tournament(std::vector<Individual>& popArray) {

	// Randomly select two different indices within the range of the population array
	int max_index = static_cast<int>(popArray.size()) - 1;
	int choice1 = std::rand() % (max_index + 1);
	int choice2 = std::rand() % (max_index + 1);

	// Ensure the two choices are different
	while (choice1 == choice2) {
		choice2 = std::rand() % (max_index + 1);
	}

	// Since the population array is sorted by fitness (lower index is fitter),
	// select the individual with the lower index
	int chosen_index = (choice1 <= choice2) ? choice1 : choice2;

	// Set the progenator to the selected individual
	progenator_ = popArray[chosen_index];
}