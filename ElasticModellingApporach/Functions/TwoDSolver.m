%% 2D : Referential position is +X

clear; clc; %close all;

% Numeric Parameters
n = 20; % Number of discrete linear segments
L = 0.05; % Total elastica length
s = linspace(0, L, n)'; % Column vector
delta_s = repmat(s(2) - s(1), n, 1); % Local, immutable link length
theta = zeros(n, 1); % Rotation angles
theta_p = zeros(n, 1); % Derivative of theta with respect to s
r = zeros(n, 2); % Position vector as a function of s

% Geometric Parameters
radius = 2e-3; % Radius [m]
A = pi * radius^2; % Cross-sectional area [m^2]
V = A; % Volume per unit length [m^2]
I = 0.25 * pi * radius^4; % Second moment of area [m^4]

% Material Properties
E = 10e5; % Young's modulus [N/m^2]
rho = 1100; % Density [kg/m^3]
g = -9.81; % Acceleration due to gravity [m/s^2]

% Magnetic Properties
mu0 = 4 * pi * 1e-7; % Air magnetic permeability [H/m or N/A^2]
remanence = 0.10; % Saturation remanent magnetic field [T]
M = remanence * V / mu0; % Magnitude of magnetic moment [A·m^2]
m = zeros(n, 2); % Referential magnetic moment vectors
m_def = zeros(n, 2); % Deformed magnetic moment vectors
m_hat = zeros(n, 2); % Derivative of m_def with respect to theta

% Initialize Magnetic Moments
m(:, 1) = M; % Magnetic moment along the x-axis
m(1:round(n/2), 1) = M; % First half segments
%m(round(n/2)+1:end, 1) = -M; % Uncomment if needed

% Magnetic Field (Constants)
Bx_const = 0.001880; 
By_const = -0.001672; 
B_interp = repmat([Bx_const, By_const], n, 1);

% Loop Parameters
loops = 20; % Maximum number of iterations
mu = 0.5; % Numerical damping factor
theta_p_old = ones(n, 1) * 1e-6;

% Precompute Constants
const1 = (rho * g * A) / (E * I);
const2 = 1 / (E * I);
delta_s1 = delta_s(1);

% Initial Position
r(1, :) = [0, 0];

% Iterative Loop
for j = 1:loops
    % Compute cos and sin of theta
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    
    % Compute m_def and m_hat without loops
    m_def(:, 1) = cos_theta .* m(:, 1) - sin_theta .* m(:, 2);
    m_def(:, 2) = sin_theta .* m(:, 1) + cos_theta .* m(:, 2);
    m_hat(:, 1) = -sin_theta .* m(:, 1) - cos_theta .* m(:, 2);
    m_hat(:, 2) = cos_theta .* m(:, 1) - sin_theta .* m(:, 2);
    
    % Map to Cartesian Coordinates
    dx = delta_s1 * cos_theta;
    dy = delta_s1 * sin_theta;
    r(2:end, 1) = cumsum(dx(1:end-1)) + r(1, 1);
    r(2:end, 2) = cumsum(dy(1:end-1)) + r(1, 2);
    
    % Compute Dot Product Efficiently
    dot_product = m_hat(:, 1) * Bx_const + m_hat(:, 2) * By_const;
    
    % Calculate Theta Prime Based on PDE Using Real Cube Root
    expr = const1 * cos_theta + (theta_p_old.^2) .* (const2 * dot_product);
    theta_p = (1 - mu) * theta_p_old + mu * cbrt(expr);
    
    % Integrate Theta Prime for Theta
    theta(2:end) = theta(1:end-1) + ((theta_p(2:end) + theta_p(1:end-1)) / 2) * delta_s1;
    
    % Convergence Check
    err = 1e-6;
    if max(abs(theta_p - theta_p_old)) < err
        fprintf('norm of error = %e\n', norm(theta_p - theta_p_old));
        fprintf('Converged at iteration j = %d\n', j);
        break;
    end
    
    if j == loops
        fprintf('norm of error = %e\n', norm(theta_p - theta_p_old));
        fprintf('Unconverged: reached maximum iterations j = %d\n', j);
    end
    
    % Update theta_p_old
    theta_p_old = theta_p;
end

% Verify Total Length of the Elastica
dr = diff(r);
total_length = sum(sqrt(sum(dr.^2, 2)));
fprintf('Check elastica total length = %.6f\n', total_length);

% Plotting
scale = 0.001; % Scale factor for quiver plots

figure(1);
plot(real(r(:, 1)), real(r(:, 2)), 'k.-');
xlabel('x position (m)');
ylabel('y position (m)');
hold on;
% Plot Magnetic Moment Vectors
quiver(real(r(:, 1)), real(r(:, 2)), real(m_def(:, 1)) * scale, real(m_def(:, 2)) * scale, 0, 'MaxHeadSize', 0.5);
% Plot Magnetic Field Vectors (Optional)
%quiver(real(r(:, 1)), real(r(:, 2)), Bx_const * scale, By_const * scale, 0, 'MaxHeadSize', 0.5, 'Color', 'b');
plot(0, 0, 'rx');
axis equal;
title('Elastica with Magnetic Moments and Field');
hold off;

% --- Function Definitions ---

function y = cbrt(x)
    y = sign(x) .* abs(x).^(1/3);
end