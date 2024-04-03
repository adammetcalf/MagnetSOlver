close all;
clear;
clc;

% Define world inputs
JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0];
TentacleMagnetisation = 12;
LinkLength = 0.01;
HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.

% Create World
world = World(LinkLength, JointAngles, TentacleMagnetisation, HomogeneousField, MultipoleActive);

% Plot the world
world = world.plotWorld(false,false,1);

Angles2 = world.getJointAngles();

FT = world.getForcesTorques();

% Make new angles
JointAngles = [90,180;0,0;0,0;0,0;0,0;0,0];

% Update world with new angles
world = world.UpDateAngles(JointAngles);

% Plot the world
world = world.plotWorld(false,false,2);

Angles2 = world.getJointAngles();

FT2 = world.getForcesTorques();