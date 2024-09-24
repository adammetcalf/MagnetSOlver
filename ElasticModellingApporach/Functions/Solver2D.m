function [xpos, ypos] = Solver2D(n, L, Bxin, Byin, magMoments, orientation)
    % Solver2D computes the position of a flexible tentacle under magnetic and gravitational forces.
    %
    % Inputs:
    %   n           - Number of nodes along the tentacle
    %   L           - Total length of the tentacle [m]
    %   Bxin        - X-component of the magnetic field [T]
    %   Byin        - Y-component of the magnetic field [T]
    %   magMoments  - Magnetic moments at each node (n x 2 matrix)
    %   orientation - Starting orientation (0 for horizontal, 1 for vertical downwards)
    %
    % Outputs:
    %   xpos        - X positions of the tentacle nodes
    %   ypos        - Y positions of the tentacle nodes

    % Numeric Parameters
    s = linspace(0, L, n)'; % Column vector
    delta_s = s(2) - s(1); % Link length
    theta = zeros(n, 1); % Rotation angle at each node
    theta_p = zeros(n, 1); % First derivative of theta
    theta_p_old = ones(n, 1) * 1e-6; % Previous iteration of theta_p
    r = zeros(n, 2); % Position vector at each node

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
    remanence = 0.05; % Saturation remanent magnetic field [T]
    M = remanence * V / mu0; % Magnitude of magnetic moment [A·m^2]
    m = M * magMoments; % Scale the provided magnetic moments by M
    m_def = zeros(n, 2); % Deformed magnetic moment vectors
    m_hat = zeros(n, 2); % Derivative of m_def with respect to theta

    % Magnetic Field Components
    Bx_const = Bxin;
    By_const = Byin;
    B_interp = repmat([Bx_const, By_const], n, 1);

    % Loop Parameters
    loops = 100; % Maximum number of iterations
    mu = 0.5; % Numerical damping parameter

    % Precompute Constants
    EI = E * I; % Bending stiffness
    const1 = (rho * g * A) / EI;
    const2 = 1 / EI;

    % Initial Position
    r(1, :) = [0, 0]; % Starting position at the origin

    % Set Initial Orientation Based on Input Parameter
    if orientation == 0
        % Horizontal orientation
        theta = zeros(n, 1); % Rotation angle at each node
    elseif orientation == 1
        % Vertically downwards orientation
        theta = (-pi / 2) * ones(n, 1); % Rotation angle at each node set to -90 degrees
    else
        error('Invalid orientation. Use 0 for horizontal, 1 for vertically downwards.');
    end

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
        dx = delta_s * cos_theta;
        dy = delta_s * sin_theta;
        r(2:end, 1) = cumsum(dx(1:end-1)) + r(1, 1);
        r(2:end, 2) = cumsum(dy(1:end-1)) + r(1, 2);

        % Compute Dot Product Efficiently
        dot_product = m_hat(:, 1) * Bx_const + m_hat(:, 2) * By_const;

        % Calculate Theta Prime Using Real Cube Root
        expr = const1 * cos_theta + (theta_p_old.^2) .* (const2 * dot_product);
        theta_p = (1 - mu) * theta_p_old + mu * cbrt(expr);

        % Integrate Theta Prime to Get Theta
        theta(2:end) = theta(1:end-1) + ((theta_p(2:end) + theta_p(1:end-1)) / 2) * delta_s;

        % Convergence Check Using Maximum Absolute Difference
        err_tol = 1e-6;
        if max(abs(theta_p - theta_p_old)) < err_tol
            break;
        end


        % Update theta_p_old for the Next Iteration
        theta_p_old = theta_p;
    end

    % Return the x and y Positions
    xpos = real(r(:, 1));
    ypos = real(r(:, 2));
end

% --- Helper Function ---
function y = cbrt(x)
    y = sign(x) .* abs(x).^(1/3);
end
