close all;
clear;
clc;

% Define world inputs
JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0];
TentacleMagnetisation = 12;
LinkLength = 0.01;
HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction

% Create World
world = World(LinkLength, JointAngles, TentacleMagnetisation, HomogeneousField);

% Plot the world
world = world.plotWorld(false,true);

% Make new angles
JointAngles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0];

% Update world with new angles
world = world.UpDateAngles(JointAngles);

% Plot the world
world = world.plotWorld(false,true);