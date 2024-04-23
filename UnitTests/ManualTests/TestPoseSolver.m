close all;
clear;
clc;

% This tests optimisation using both the inbuilt MATLAB Genetic algorithm
% and the inbuilt minimisation (optimisation) functions

% Define world inputs
JointAngles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
Magdirection = [0,0,0,0,0,0,0;0,0,0,0,0,0,0;1,1,1,-1,-1,-1,-1];
TentacleMagnetisation = 14;
LinkLength = 0.01;

% Spawn tentacle
tentacle = Tentacle(LinkLength,JointAngles,TentacleMagnetisation,Magdirection);
HomogeneousField = [25e-3,0,0]; % Magnetic field strength in Tesla (25 mT) in +z direction
MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.

% Create World
World = World(tentacle,HomogeneousField,[], MultipoleActive);

% Plot the world
World = World.plotWorld(false,false,1);

FT = World.getForcesTorques();

% Create PoseSolver
%PoseSolver = PoseSolver(World);

% Create GA 
GA = GAMatlabFunction(World);

optimizedAngles = GA.optimizeJoints();

World = World.UpDateAngles(optimizedAngles);

FT2 = World.getForcesTorques();

% Plot the world
World = World.plotWorld(false,false,2);

% Create PoseSolver
PoseSolver = PoseSolver(World);

optimizedAngles = PoseSolver.optimizeJoints();

World = World.UpDateAngles(optimizedAngles);

% Plot the world
World = World.plotWorld(false,false,3);
