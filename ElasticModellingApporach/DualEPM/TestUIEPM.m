close all;
clear;
clc;

% Parameters
n = 20; % Number of nodes
L = 0.05; % Length of the tentacle

% Get the target positions from the UI (deforms only in x-z)
[x_positionsUI, z_positionsUI, magn_momentsUI] = elastica_UI(n, L);

% No change in y - positions allowed by the solver model.
y_positionsUI = zeros(n,1);

% Define meshgrid over which to compute and visualize the magnetic field
[x_grid, y_grid, z_grid] = meshgrid(linspace(-0.1, 0.1, 25), ...
                                    linspace(-0.1, 0.1, 25), ...
                                    linspace(-0.1, 0.1, 25));

% Define robots
robot1 = RobotWithMagnet('urdf/kuka_iiwa_1.urdf');
robot2 = RobotWithMagnet('urdf/kuka_iiwa_2.urdf');

% Robots cell array
robots = {robot1, robot2};

% Initial joint angles (concatenate the joint angles of both robots)
initialJointAngles = [robot1.getCurrentConfig(), robot2.getCurrentConfig()];

% Joint limits (concatenate the joint limits of both robots)
jointLimits1 = deg2rad(robot1.getJointLimits());
jointLimits2 = deg2rad(robot2.getJointLimits());
lb = [jointLimits1(:,1); jointLimits2(:,1)];
ub = [jointLimits1(:,2); jointLimits2(:,2)];

% Optimization options
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
                       'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 10000);

% Optimization
[optimalJointAngles, fval] = fmincon(@(jointAngles) objective_function_robot(...
    jointAngles, robots, x_positionsUI, z_positionsUI, n, L, magn_momentsUI), ...
    initialJointAngles, [], [], [], [], lb, ub, [], options);

% Update robots with the optimal joint angles
n_joints = length(optimalJointAngles) / 2;
optimalJointAngles1 = optimalJointAngles(1:n_joints);
optimalJointAngles2 = optimalJointAngles(n_joints+1:end);
robots{1} = robots{1}.updateConfiguration(optimalJointAngles1);
robots{2} = robots{2}.updateConfiguration(optimalJointAngles2);

% Visualization of robots and tentacle
figure;
hold on;
robots{1}.showRobot();
robots{2}.showRobot();

% Define the magnetic field function handle for the tentacle
magFieldFunc = @(x_positions, z_positions) computeCombinedMagneticFieldAtTentacle(robots, x_positions, z_positions);

% Solve the elastica using the Solver2DTwo function with the computed magnetic field
[x_positions, z_positions] = Solver2DTwo(n, L, magFieldFunc, magn_momentsUI);

% Plot the resulting tentacle shape
plot3(x_positions, y_positionsUI ,z_positions, 'k.-', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Optimized Tentacle Shape with Robots');
axis equal;
hold off;

% Visualize the magnetic field over the meshgrid
[Bx_total, By_total, Bz_total] = computeCombinedMagneticFieldAtGrid(robots, x_grid, y_grid, z_grid);

% Plot the magnetic field
figure;
quiver3(x_grid, y_grid, z_grid, Bx_total, By_total, Bz_total);
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Magnetic Field Visualization');
axis equal;

%% --- Functions ---

function error = objective_function_robot(jointAngles, robots, x_positionsUI, z_positionsUI, n, L, magMoments)
    % Extract joint angles for each robot
    n_joints = length(jointAngles) / 2;
    jointAngles1 = jointAngles(1:n_joints);
    jointAngles2 = jointAngles(n_joints+1:end);

    % Update robots' configurations
    robots{1} = robots{1}.updateConfiguration(jointAngles1);
    robots{2} = robots{2}.updateConfiguration(jointAngles2);

    % Define the magnetic field function handle
    magFieldFunc = @(x_positions, z_positions) computeCombinedMagneticFieldAtTentacle(robots, x_positions, z_positions);

    % Solve the elastica using the Solver2DTwo function with the computed magnetic field
    [x_positions, z_positions] = Solver2DTwo(n, L, magFieldFunc, magMoments);

    % Calculate the error (sum of squared distances)
    error = sum((x_positions - x_positionsUI).^2 + (z_positions - z_positionsUI).^2);
end

function [Bx_total, Bz_total] = computeCombinedMagneticFieldAtTentacle(robots, x_positions, z_positions)
    Bx_total = zeros(size(x_positions));
    Bz_total = zeros(size(z_positions));

    % Loop over each robot
    for k = 1:length(robots)
        [Bx, Bz] = computeMagneticFieldAtPositions(robots{k}, x_positions, zeros(size(x_positions)), z_positions);
        Bx_total = Bx_total + Bx;
        Bz_total = Bz_total + Bz;
    end
end

function [Bx_total, By_total, Bz_total] = computeCombinedMagneticFieldAtGrid(robots, x_grid, y_grid, z_grid)
    Bx_total = zeros(size(x_grid));
    By_total = zeros(size(y_grid));
    Bz_total = zeros(size(z_grid));

    % Loop over each robot
    for k = 1:length(robots)
        [Bx, By, Bz] = computeMagneticFieldAtPositions(robots{k}, x_grid(:), y_grid(:), z_grid(:));
        % Reshape to grid size
        Bx = reshape(Bx, size(x_grid));
        By = reshape(By, size(y_grid));
        Bz = reshape(Bz, size(z_grid));
        Bx_total = Bx_total + Bx;
        By_total = By_total + By;
        Bz_total = Bz_total + Bz;
    end
end

function [Bx, By, Bz] = computeMagneticFieldAtPositions(robot, x_positions, y_positions, z_positions)
    % Get the magnet's position and magnetic moment
    [magnetPosition, magnetMoment] = robot.getMagneticField();

    % Positions where we compute the field
    positions = [x_positions(:), y_positions(:), z_positions(:)]; % Nx3 matrix

    % Initialize magnetic field arrays
    Bx = zeros(size(x_positions));
    By = zeros(size(y_positions));
    Bz = zeros(size(z_positions));

    % Magnetic field constants
    mu0 = 4 * pi * 1e-7; % Vacuum permeability

    % Compute the magnetic field at each point
    for i = 1:length(x_positions)
        % Position vector from magnet to field point
        r_vec = positions(i,:)' - magnetPosition; % 3x1 vector
        r = norm(r_vec);
        if r == 0
            continue; % Avoid singularity
        end
        r_hat = r_vec / r;

        % Compute the magnetic field using the magnetic dipole formula
        m = magnetMoment;
        B = (mu0 / (4 * pi * r^3)) * (3 * (m' * r_hat) * r_hat - m);

        Bx(i) = B(1);
        By(i) = B(2);
        Bz(i) = B(3);
    end
end