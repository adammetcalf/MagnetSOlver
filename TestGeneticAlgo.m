close all;
clear;
clc;

% Define world inputs
JointAngles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0];
Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
TentacleMagnetisation = 12;
LinkLength = 0.01;
HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.


% Create World
world = World(LinkLength, JointAngles, TentacleMagnetisation, Magdirection, HomogeneousField,[], MultipoleActive);

% Plot the world
world = world.plotWorld(false,false,1);

GeneticAlgo = GeneticAlgorithmSolution(world);

BestAngles = GeneticAlgo.Train();

% Update world with new angles
world = world.UpDateAngles(BestAngles);

% Plot the world
world = world.plotWorld(false,false,2);

FT = world.getForcesTorques();

