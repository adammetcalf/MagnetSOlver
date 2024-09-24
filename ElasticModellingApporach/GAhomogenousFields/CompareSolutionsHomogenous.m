close all;
clear;
clc;

% Parameters
n = 50; % Number of nodes
L = 0.04; % Length of the tentacle
orientation =1;

% Get the target positions from the UI
[x_positionsUI, y_positionsUI, magn_momentsUI] = elastica_UI(n, L);

%% MATLAB Optimisation Algorithm (fmincon)
% Initial guess for the magnetic field [Bx, By]
initial_field = [0.001, 0.001];

% Set lower and upper bounds for Bx and By (adjust as needed)
lb = [-0.0025, -0.0025]; % Lower bounds
ub = [0.0025, 0.0025];   % Upper bounds

% Set optimization options
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
                       'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 500);

% Benchmark start
tic;
% Use fmincon with the objective function
optimal_field = fmincon(@(field) objective_function(field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation), ...
                        initial_field, [], [], [], [], lb, ub, [], options);
% Benchmark end
fmincon_time = toc;

% Solve with the optimal field
[xpos, ypos] = Solver2D(n, L, optimal_field(1), optimal_field(2), magn_momentsUI,orientation);

% Calculate the error for fmincon
fmincon_error = objective_function(optimal_field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation);

% Output the optimal field
fprintf('fmincon - Optimal Bx: %.6f, By: %.6f\n', optimal_field(1), optimal_field(2));
fprintf('fmincon - Execution Time: %.6f seconds, Error: %.6f\n', fmincon_time, fmincon_error);

%% Genetic Algorithm
% Create world
world = world(x_positionsUI, y_positionsUI, magn_momentsUI,L,orientation);

% Instantiate GA
GAsolution = GA(world);

% Benchmark start
tic;
% Train GA
optimal_field2 = GAsolution.Train(true);
% Benchmark end
ga_time = toc;

% Solve with the optimal field
[xpos2, ypos2] = Solver2D(n, L, optimal_field2(1), optimal_field2(2), magn_momentsUI,orientation);

% Calculate the error for GA
ga_error = objective_function(optimal_field2, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation);

% Output the optimal field
fprintf('GA - Optimal Bx: %.6f, By: %.6f\n', optimal_field2(1), optimal_field2(2));
fprintf('GA - Execution Time: %.6f seconds, Error: %.6f\n', ga_time, ga_error);

%% Particle Swarm Optimization (PSO)
% Define the objective function
objectiveFunction = @(field) objective_function(field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation);

% Set options for particleswarm
options = optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 100, 'Display', 'iter');

% Benchmark start
tic;
% Run PSO
[optimal_field3, fval] = particleswarm(objectiveFunction, 2, lb, ub, options);
% Benchmark end
pso_time = toc;

% Solve with the optimal field
[xpos3, ypos3] = Solver2D(n, L, optimal_field3(1), optimal_field3(2), magn_momentsUI,orientation);

% Calculate the error for PSO
pso_error = objective_function(optimal_field3, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation);

% Output the optimal field
fprintf('PSO - Optimal Bx: %.6f, By: %.6f\n', optimal_field3(1), optimal_field3(2));
fprintf('PSO - Execution Time: %.6f seconds, Error: %.6f\n', pso_time, pso_error);


%% MATLAB GA

% Set lower and upper bounds for Bx and By (adjust as needed)
lb = [-0.0025, -0.0025]; % Lower bounds
ub = [0.0025, 0.0025];   % Upper bounds

%Define the objective function
objectiveFunction = @(field) objective_function(field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation);

% Genetic Algorithm options setup
options = optimoptions('ga', ...
                               'PopulationSize', 100, ...
                               'MaxGenerations', 100, ...
                               'CrossoverFcn', @crossoverintermediate, ...
                               'MutationFcn', @mutationadaptfeasible, ...
                               'Display', 'iter');  % Iterative display of the process

%Benchmark start
tic;
% Use ga with the objective function
optimal_field_ga = ga(objectiveFunction, 2, [], [], [], [], lb, ub, [], options);
% Benchmark end
MATLAB_GA_time = toc;

% Solve with the optimal field
[xpos_ga, ypos_ga] = Solver2D(n, L, optimal_field_ga(1), optimal_field_ga(2), magn_momentsUI,orientation);

% Calculate the error for MATLAB GA
matlab_ga_error = objective_function(optimal_field_ga, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation);

% Output the optimal field
fprintf('MATLAB GA - Optimal Bx: %.6f, By: %.6f\n', optimal_field_ga(1), optimal_field_ga(2));
fprintf('MATLAB GA - Execution Time: %.6f seconds, Error: %.6f\n', MATLAB_GA_time, matlab_ga_error);

%% Plotting Results

% Plotting the desired shape
figure(1);
plot(x_positionsUI, y_positionsUI, 'k.-');
hold on;
plot(0, 0, 'rx');
hold off;
xlabel('x position (m)');
ylabel('y position (m)');
axis equal;
title('Desired Elastica Shape');

% Plotting the optimized shape (fmincon)
figure(2);
plot(xpos, ypos, 'k.-');
hold on;
plot(0, 0, 'rx');
hold off;
xlabel('x position (m)');
ylabel('y position (m)');
axis equal;
title('Elastica Shape with fmincon Field');

% Plotting the optimized shape (GA)
figure(3);
plot(xpos2, ypos2, 'k.-');
hold on;
plot(0, 0, 'rx');
hold off;
xlabel('x position (m)');
ylabel('y position (m)');
axis equal;
title('Elastica Shape with GA Field');

% Plotting the optimized shape (PSO)
figure(4);
plot(xpos3, ypos3, 'k.-');
hold on;
plot(0, 0, 'rx');
hold off;
xlabel('x position (m)');
ylabel('y position (m)');
axis equal;
title('Elastica Shape with PSO Field');

% Plotting the optimized shape (MATLAB GA)
figure(5);
plot(xpos_ga, ypos_ga, 'k.-');
hold on;
plot(0, 0, 'rx');
hold off;
xlabel('x position (m)');
ylabel('y position (m)');
axis equal;
title('Elastica Shape with MATLAB GAField');


%% Objective Function
function error = objective_function(field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation)
    % Set the magnetic field
    Bx = field(1);
    By = field(2);

    % Solve the elastica using the Solver2D function
    [x_positions, y_positions] = Solver2D(n, L, Bx, By, magn_momentsUI,orientation);

    % Calculate the error (Euclidean distance)
    error = sqrt(sum((x_positions - x_positionsUI).^2 + (y_positions - y_positionsUI).^2));
end

