%% 2D : Referential position is +X

clear; clc; close all;


% Numeric
n = 20; % Number of discrete linear segments
L = 0.05; % Total elastica length
s = linspace(0, L, n)'; % Column vector
delta_s = repmat(s(2) - s(1), n, 1); % Local, immutable, Link Length
theta = zeros(n, 1); % Rotation Matrix from d3 to e3
theta_p = zeros(n, 1); % Differential of R wrt s
r = zeros(n, 2); % Vector position as a function of "s"

% Geometric
radius = 2e-3; % Radius [m]
A = pi * radius^2; % Cross Sectional Area [m^2]
V = A; % Volume per unit length [m^2]
I = 0.25 * pi * radius^4; % Second moment of area [m^4]

% Material
E = 10e5; % Young's modulus [N/m^2]
rho = 1100; % Density [kg/m^3]
g = -9.81; % Acceleration due to gravity [m/s^2]

% Magnetic
mu0 = 4 * pi * 1e-7; % Air magnetic permeability (H/m or N/A^2)
remanence = 0.10; % Saturation Remanent Magnetic Field [T]
M = remanence * V / mu0; % Magnitude of magnetic moment (A·m^2)
m = zeros(n, 2); % Spatial vector of referential magnetic moment
m_def = zeros(n, 2); % Spatial vector of deformed magnetic moment
m_hat = zeros(n, 2); % Differential wrt theta of deformed magnetic moment

% Initialize magnetic moments
m(:, 1) = M; % Initial magnetic moment along x-axis
m(1:round(n/2), 1) = -M; % Flip direction for first half
%m(round(n/2)+1:end, 1) = M; % Second half remains positive
m(1:round(n/2), 1) = M; % NO Flip direction for first half
%m(round(n/2)+1:end, 1) = -M; % Flip direction for second half

% Create the magnetic field template
xLow = -0.1;   yLow = -0.1;
xHigh = 0.1;   yHigh = 0.1;
xIncr = 0.005;  yIncr = 0.005; % Step size

% Create the meshgrid used for magnetic simulation
[Fieldx, Fieldy] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh);

% Define the magnetic field Bx and By over the grid
%Bx = 0.001880 * ones(size(Fieldx)); 
%By = -0.001672 * ones(size(Fieldy)); 
Bx = 0.0 * ones(size(Fieldx)); 
By = -0.0 * ones(size(Fieldy));

% Loop parameters
loops = 1000; % Maximum number of iterations
mu = 0.5;
theta_p_old = ones(n, 1) * 1e-6;

R = zeros(2, 2, n); % Rotation matrices at each datapoint
R_hat = zeros(2, 2, n); % Differential of Rotation matrix wrt theta

r(1,:) = [0, 0]; % Starting position

for j = 1:loops
    % Compute rotation matrices and their derivatives
    for i = 1:n
        % Rotation Matrix for current theta
        R(:, :, i) = [cos(theta(i)), -sin(theta(i));
                      sin(theta(i)),  cos(theta(i))];
        m_def(i, :) = (R(:, :, i) * m(i, :)')';

        % Differential of Rotation Matrix wrt theta
        R_hat(:, :, i) = [-sin(theta(i)), -cos(theta(i));
                           cos(theta(i)), -sin(theta(i))];
        m_hat(i, :) = (R_hat(:, :, i) * m(i, :)')';
    end

    % Map to Cartesian coordinates
    dx = delta_s(1) * cos(theta);
    dy = delta_s(1) * sin(theta);
    r(1, :) = [0, 0];
    r(2:end, 1) = cumsum(dx(1:end-1)) + r(1,1);
    r(2:end, 2) = cumsum(dy(1:end-1)) + r(1,2);

    % Interpolate Bx and By at positions r(:,1), r(:,2)
    Bx_interp = interp2(Fieldx, Fieldy, Bx, real(r(:,1)), real(r(:,2)), 'linear', 0);
    By_interp = interp2(Fieldx, Fieldy, By, real(r(:,1)), real(r(:,2)), 'linear', 0);
    B_interp = [Bx_interp, By_interp];

    % Calculate Theta Prime based on PDE using real cube root
    expr = (rho * g * A / (E * I)) * cos(theta) + (theta_p_old.^2) .* (sum(m_hat .* B_interp, 2) / (E * I));
    theta_p = (1 - mu) * theta_p_old + mu * cbrt(expr);

    % Integrate Theta Prime for Theta
    theta(2:end) = theta(1:end-1) + ((theta_p(2:end) + theta_p(1:end-1)) / 2) * delta_s(1);

    % Convergence Clause
    err = 1e-6;
    if norm(theta_p - theta_p_old) < err
        fprintf('norm of error = %e\n', norm(theta_p - theta_p_old));
        fprintf('Converged at iteration j = %d\n', j);
        break;
    end

    if j == loops
        fprintf('norm of error = %e\n', norm(theta_p - theta_p_old));
        fprintf('Unconverged: reached maximum iterations j = %d\n', j);
    end

    % Update theta_p_old for convergence check and numerical damping
    theta_p_old = theta_p;
end

% Check elastica total length
dr = diff(r);
total_length = sum(sqrt(sum(dr.^2, 2)));
fprintf('Check elastica total length = %.6f\n', total_length);

% Plotting
figure(1);
plot(real(r(:, 1)), real(r(:, 2)), 'k.-');
hold on
plot(0,0,'rx');
xlabel('x position (m)');
ylabel('y position (m)');
axis equal;
title('Elastica Shape');
hold off;

% Adjust the scale for the quiver plot if necessary
scale = 0.001; % Adjust scale factor as needed

figure(2);
plot(real(r(:, 1)), real(r(:, 2)), 'k.-');
xlabel('x position (m)');
ylabel('y position (m)');
hold on;
% Plot the magnetic moment vectors at the correct positions
quiver(real(r(:, 1)), real(r(:, 2)), real(m_def(:, 1))*scale, real(m_def(:, 2))*scale, 0, 'MaxHeadSize', 0.5);
% Plot the magnetic field vectors at the elastica positions (optional)
quiver(real(r(:,1)), real(r(:,2)), Bx_interp*scale, By_interp*scale, 0, 'MaxHeadSize', 0.5, 'Color', 'b');
plot(0,0,'rx');
axis equal;
title('Elastica with Magnetic Moments and Field');
hold off;

% --- Function Definitions ---

function y = cbrt(x)
    y = sign(x) .* abs(x).^(1/3);
end
