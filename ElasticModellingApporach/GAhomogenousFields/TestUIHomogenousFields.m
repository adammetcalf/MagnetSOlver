close all;
clear;
clc;

% Parameters
n = 50; % Number of nodes
L = 0.05; % Length of the tentacle
Opt = 3;  % 1 = optimisation, 2 = Genetic Algorithm, 3 = particle swarm

% Get the target positions from the UI
[x_positionsUI, y_positionsUI, magn_momentsUI] = elastica_UI(n, L);

switch Opt

    case 1 % MATLAB optimisation Algorithm

    % Initial guess for the magnetic field [Bx, By]
    initial_field = [0.001, 0.001];
    
    % Set lower and upper bounds for Bx and By (adjust as needed)
    lb = [-0.0025, -0.0025]; % Lower bounds (no bounds in this case)
    ub = [0.0025, 0.0025];   % Upper bounds (no bounds in this case)
    
    % Set optimization options
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
                           'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 500);
    
    % Use fmincon with the modified objective function
    optimal_field = fmincon(@(field) objective_function(field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI), ...
                            initial_field, [], [], [], [], lb, ub, [], options);
    
    % Output the optimal field
    fprintf('Optimal Bx: %.6f, By: %.6f\n', optimal_field(1), optimal_field(2));
    
    % Solve with the optimal field
    [xpos, ypos] = Solver2D(n, L, optimal_field(1), optimal_field(2), magn_momentsUI);

    case 2 % Genetic Algorithm

    % Create world
    world = world(x_positionsUI, y_positionsUI, magn_momentsUI,L);

    % Instantiate GA
    GAsolution = GA(world);

    % Train GA
    optimal_field = GAsolution.Train(true);

    % Solve with the optimal field
    [xpos, ypos] = Solver2D(n, L, optimal_field(1), optimal_field(2), magn_momentsUI);

    case 3 % particle swarm
    
    % Define the objective function
    objectiveFunction = @(field) objective_function(field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI);
    
    % Set lower and upper bounds for [Bx, By]
    lb = [-0.0025, -0.0025];
    ub = [0.0025, 0.0025];
    
    % Set options for particleswarm
    options = optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 100, 'Display', 'iter');
    
    % Run PSO
    [optimal_field, fval] = particleswarm(objectiveFunction, 2, lb, ub, options);

    % Solve with the optimal field
    [xpos, ypos] = Solver2D(n, L, optimal_field(1), optimal_field(2), magn_momentsUI);

end



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

% Plotting the optimized shape
figure(2);
plot(xpos, ypos, 'k.-');
hold on;
plot(0, 0, 'rx');
hold off;
xlabel('x position (m)');
ylabel('y position (m)');
axis equal;
title('Elastica Shape with Optimized Field');

%% Functions

function error = objective_function(field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI)
    % Set the magnetic field
    Bx = field(1);
    By = field(2);

    % Solve the elastica using the Solver2D function
    [x_positions, y_positions] = Solver2D(n, L, Bx, By, magn_momentsUI);

    % Calculate the error (Euclidean distance)
    error = sqrt(sum((x_positions - x_positionsUI).^2 + (y_positions - y_positionsUI).^2));
end
