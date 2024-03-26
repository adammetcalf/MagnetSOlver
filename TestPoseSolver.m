close all;
clear;
clc;

% Define world inputs
JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
TentacleMagnetisation = 12;
LinkLength = 0.01;
HomogeneousField = [0,0,-25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.

% Create World
World = World(LinkLength, JointAngles, TentacleMagnetisation, HomogeneousField, MultipoleActive);

% Plot the world
World = World.plotWorld(false,false);

FT = World.getForcesTorques();

% Create PoseSolver
PoseSolver = PoseSolver(World);

optimizedAngles = PoseSolver.optimizeJoints();

World = World.UpDateAngles(optimizedAngles);

FT2 = World.getForcesTorques();

% Plot the world
World = World.plotWorld(false,false);


