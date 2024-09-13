close all;
clear;
clc;

%%  Define Tentacle using UI
% Declare the global variable to retrieve the result
global tentacleResult;
global BestAngles;

% Initialize or clear the global variable
tentacleResult = [];

% open the GUI
fig = simpleApp; 

% Wait for the UI to close
waitfor(fig);

%%  define world inputs
HomogeneousField = [25e-3,0,0]; % Magnetic field strength in Tesla (25 mT) in +z direction
MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.

%%  Define EPM
location = [0.15;0;0];
location2 = [0;0;-0.15];
orientation = [1,0,0;0,1,0;0,0,1];
magnet = ExternalEPM(location,orientation);
magnet2 = ExternalEPM(location2,orientation);


% Create World
world = World(tentacleResult, HomogeneousField,[], MultipoleActive);


% Make new angles
JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0];

% Update world with new angles
world = world.UpDateAngles(JointAngles);

% Plot the world
world = world.plotWorld(false,false,1);

FT = world.getForcesTorques();

% Make new angles
JointAngles = [90,270;0,0;0,0;0,0;0,0;0,0];

% Update world with new angles
world = world.UpDateAngles(JointAngles);

% Plot the world
world = world.plotWorld(false,false,2);

FT2 = world.getForcesTorques();