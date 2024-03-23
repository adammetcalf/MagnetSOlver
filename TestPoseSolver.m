close all;
clear;
clc;

% Define world inputs
JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0];
TentacleMagnetisation = 12;
LinkLength = 0.01;
HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
MultipoleActive = true; % Include the magnetic effects of the tentacle links on the magentic field.

% Create World
world = World(LinkLength, JointAngles, TentacleMagnetisation, HomogeneousField, MultipoleActive);

% Create PoseSolver
PoseSolver = PoseSolver(world);

optimizedAngles = PoseSolver.optimizeJoints()
