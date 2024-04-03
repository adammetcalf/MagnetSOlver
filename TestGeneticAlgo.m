close all;
clear;
clc;

% Define world inputs
JointAngles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0];
TentacleMagnetisation = 12;
LinkLength = 0.01;
HomogeneousField = [25e-3,0,0]; % Magnetic field strength in Tesla (25 mT) in +z direction
MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.

% Create World
world = World(LinkLength, JointAngles, TentacleMagnetisation, HomogeneousField, MultipoleActive);

% Plot the world
world = world.plotWorld(false,false,1);

GeneticAlgo = GeneticAlgorithmSolution(world);

[GeneticAlgo, BestAngles] = GeneticAlgo.Train();

% Update world with new angles
world = world.UpDateAngles(BestAngles);

% Plot the world
world = world.plotWorld(false,false,2);

FT = world.getForcesTorques();

% Create PoseSolver
PoseSolver = PoseSolver(world);

optimizedAngles = PoseSolver.optimizeJoints();

world = world.UpDateAngles(optimizedAngles);

% Plot the world
world = world.plotWorld(false,false,3);

FT2 = world.getForcesTorques();
