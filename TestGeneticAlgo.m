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


% Plot the world
world = world.plotWorld(false,false,1);

%%%%% Lgecay - from before UI integration
%GeneticAlgo = GeneticAlgorithmSolution(world);
%BestAngles = GeneticAlgo.Train();

fig = GAapp(world);

% Wait for the UI to close
waitfor(fig);

% Update world with new angles
world = world.UpDateAngles(BestAngles);

% Plot the world
world = world.plotWorld(false,false,2);

FT = world.getForcesTorques();


% #TODO - Stiffness is really fucking up the solution

% #TODO - something about the way we are dealing with the individuals after
% optimisation seems to be having the same effect on the solution as if
% gravity is upside down. oSomething to do with inverting the best
% individudal??

% #TODO optimise theta angles too
% #TODO integrate the alpha2 angles into the solution


% #TODO tentacle magentisation should be defined by the mass of particle
% inclusions

% #TODO angle constraints in optimisation and genetic (individual) must be
% reconsidered (particulary 3d)

% #TODO The iron should be considered - magnetisation/ rememnace?


% It makes sense that we seem to get results that seem dierectly opposite
% to what we expect because of the local minim resilt where we have no
% torques to rtate us into the expected result direction

% #TODO IK
