% main.m

close all;
clear;
clc;

% Parameters
n = 20; % Number of nodes
L = 0.05; % Length of the tentacle

% Get the target positions from the UI
[x_positionsUI, z_positionsUI, magn_momentsUI] = elastica_UI(n, L);


% Define robots
robot1_template = RobotWithMagnet('urdf/kuka_iiwa_1.urdf');
robot2_template = RobotWithMagnet('urdf/kuka_iiwa_2.urdf');

% Get initial configurations
robot1_currentConfig = robot1_template.getCurrentConfig();
robot2_currentConfig = robot2_template.getCurrentConfig();

% Get joint limits
robot1_jointLimits = robot1_template.getJointLimits(); % size (7 x 2)
robot2_jointLimits = robot2_template.getJointLimits();

% Flatten joint limits
lb_robot1 = robot1_jointLimits(:,1)';
ub_robot1 = robot1_jointLimits(:,2)';
lb_robot2 = robot2_jointLimits(:,1)';
ub_robot2 = robot2_jointLimits(:,2)';

% Combine initial configurations and joint limits
initial_joint_angles = [robot1_currentConfig, robot2_currentConfig];
lb = [lb_robot1, lb_robot2];
ub = [ub_robot1, ub_robot2];

% Set optimization options
%options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
%                       'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 1000);

% Set options for particleswarm
options = optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 100, 'Display', 'iter', 'UseParallel', false);

% Define the objective function as a function handle
objectiveFunction = @(jointAngles) objective_function(jointAngles, x_positionsUI, z_positionsUI, n, L, magn_momentsUI, robot1_template, robot2_template);

% Number of variables
numVars = length(initial_joint_angles); % Should be 14

% Run fmincon
%[optimal_joint_angles, fval] = fmincon(objectiveFunction, initial_joint_angles, [], [], [], [], lb, ub, [], options);

% Run particleswarm
[optimal_joint_angles, fval] = particleswarm(objectiveFunction, numVars, lb, ub, options);

% Update robots with optimal joint angles
optimal_jointAngles_robot1 = optimal_joint_angles(1:7);
optimal_jointAngles_robot2 = optimal_joint_angles(8:14);

robot1 = robot1_template.updateConfiguration(optimal_jointAngles_robot1);
robot2 = robot2_template.updateConfiguration(optimal_jointAngles_robot2);

% Compute the final magnetic field function
magFieldFunc = @(x_positions, z_positions) compute_magnetic_field(x_positions, z_positions, robot1, robot2);

% Solve with the optimal field
[xpos, zpos] = Solver2DTwo(n, L, magFieldFunc, magn_momentsUI);

% Plotting the desired shape
figure(1);
plot(x_positionsUI, z_positionsUI, 'k.-');
hold on;
plot(0, 0, 'rx');
hold off;
xlabel('x position (m)');
ylabel('z position (m)');
axis equal;
title('Desired Tentacle Shape');

% Plotting the optimized shape
figure(2);
plot(xpos, zpos, 'k.-');
hold on;
plot(0, 0, 'rx');
hold off;
xlabel('x position (m)');
ylabel('z position (m)');
axis equal;
title('Tentacle Shape with Optimized Robot Configurations');

% Visualization of robots and magnetic field
figure(3);
hold on;
robot1.showRobot();
robot2.showRobot();
plot3(xpos, zeros(n,1), zpos, 'k.-', 'LineWidth', 2);
xlabel('x');
ylabel('y');
zlabel('z');
title('Robots with Optimal Configurations and Magnetic Field');
view(3);
grid on;

% --- Plot the magnetic field onto figure 3 ---
% Define grid over x and z
x_min = min(xpos) - 0.01;
x_max = max(xpos) + 0.01;
z_min = min(zpos) - 0.01;
z_max = max(zpos) + 0.01;

num_points_x = 20;
num_points_z = 20;

x_range = linspace(x_min, x_max, num_points_x);
z_range = linspace(z_min, z_max, num_points_z);

[X_grid, Z_grid] = meshgrid(x_range, z_range);

% Flatten the grids
X_flat = X_grid(:);
Z_flat = Z_grid(:);

% Initialize B components
Bx_flat = zeros(size(X_flat));
By_flat = zeros(size(X_flat));
Bz_flat = zeros(size(X_flat));

% Get the magnet positions and moments from the robots
[position1, moment1] = robot1.getMagneticField();
[position2, moment2] = robot2.getMagneticField();

% Compute the magnetic field at each point
for i = 1:length(X_flat)
    point = [X_flat(i); 0; Z_flat(i)]; % y = 0

    % Magnetic field from robot1
    B1 = calculateMagneticField(point, position1, moment1);

    % Magnetic field from robot2
    B2 = calculateMagneticField(point, position2, moment2);

    % Total magnetic field
    B_total = B1 + B2;

    % Store B components
    Bx_flat(i) = B_total(1);
    By_flat(i) = B_total(2);
    Bz_flat(i) = B_total(3);
end

% Reshape B components back to grid shape
Bx_grid = reshape(Bx_flat, size(X_grid));
By_grid = reshape(By_flat, size(X_grid));
Bz_grid = reshape(Bz_flat, size(X_grid));

% Plot the magnetic field vectors
% For visualization, scale the vectors appropriately
scale = 1e6; % Adjust scale for better visualization
quiver3(X_grid, zeros(size(X_grid)), Z_grid, Bx_grid*scale, By_grid*scale, Bz_grid*scale, 'r');

hold off;

%% Objective Function
function error = objective_function(jointAngles, x_positionsUI, z_positionsUI, n, L, magn_momentsUI, robot1_template, robot2_template)

    % Create copies of the robot objects to avoid modifying shared state
    robot1 = robot1_template;
    robot2 = robot2_template;

    % Extract joint angles for each robot
    jointAngles_robot1 = jointAngles(1:7);
    jointAngles_robot2 = jointAngles(8:14);

    % Update robot configurations
    robot1 = robot1.updateConfiguration(jointAngles_robot1);
    robot2 = robot2.updateConfiguration(jointAngles_robot2);

    % Define the magnetic field function
    magFieldFunc = @(x_positions, z_positions) compute_magnetic_field(x_positions, z_positions, robot1, robot2);

    % Solve the elastica using Solver2DTwo
    [x_positions, z_positions] = Solver2DTwo(n, L, magFieldFunc, magn_momentsUI);

    % Calculate the error (Euclidean distance)
    error = sqrt(sum((x_positions - x_positionsUI).^2 + (z_positions - z_positionsUI).^2));

end

%% Compute Magnetic Field Function
function [Bx_interp, Bz_interp] = compute_magnetic_field(x_positions, z_positions, robot1, robot2)

    num_points = length(x_positions);
    Bx = zeros(num_points, 1);
    Bz = zeros(num_points, 1);

    % Get the magnet positions and moments from the robots
    [position1, moment1] = robot1.getMagneticField();
    [position2, moment2] = robot2.getMagneticField();

    for i = 1:num_points
        point = [x_positions(i); 0; z_positions(i)]; % y = 0

        % Magnetic field from robot1
        B1 = calculateMagneticField(point, position1, moment1);

        % Magnetic field from robot2
        B2 = calculateMagneticField(point, position2, moment2);

        % Total magnetic field
        B_total = B1 + B2;

        % Store Bx and Bz
        Bx(i) = B_total(1);
        Bz(i) = B_total(3); % Since z is the third coordinate
    end

    % Return Bx_interp and Bz_interp
    Bx_interp = Bx;
    Bz_interp = Bz;

end

%% Calculate Magnetic Field Function
function B_total = calculateMagneticField(point, position, moment)

    % Compute displacement vector R from dipole to point
    R = point - position; % Both are column vectors

    % Compute the norm of R
    R_norm = norm(R);

    % Compute dot product of moment and R
    m_dot_R = dot(moment, R);

    % Compute the magnetic field using dipole formula
    if R_norm ~= 0
        B_total = (3 * R * m_dot_R / R_norm^5) - (moment / R_norm^3);
    else
        B_total = [0;0;0]; % Handle singularity
    end

end
