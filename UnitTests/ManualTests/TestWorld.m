close all;
clear;
clc;

% Define world inputs
JointAngles = [90,180;0,0;0,0;0,0;0,0;0,0];
Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
TentacleMagnetisation = 12;
LinkLength = 0.01;
HomogeneousField = [0,0,0]; % Magnetic field strength in Tesla (25 mT) in +z direction
MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.

% Define EPM
location = [0.15;0;0];
location2 = [-0.15;0;0];
orientation = [1,0,0;0,1,0;0,0,1];
magnet = ExternalEPM(location,orientation);
magnet2 = ExternalEPM(location2,orientation);

% Create World
world = World(LinkLength, JointAngles, TentacleMagnetisation, Magdirection, HomogeneousField, [magnet,magnet2], MultipoleActive);

% Plot the world
world = world.plotWorld(false,true,1);

Angles2 = world.getJointAngles();

FT = world.getForcesTorques();

% Make new angles
JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0];

% Update world with new angles
world = world.UpDateAngles(JointAngles);

% Plot the world
world = world.plotWorld(false,true,2);

Angles2 = world.getJointAngles();

FT2 = world.getForcesTorques();