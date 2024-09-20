function [xpos, ypos] = Solver2D(n, L, Bxin, Byin, magMoments)
    % Numeric
    s = linspace(0, L, n)'; % Column vector
    delta_s = repmat(s(2) - s(1), n, 1); % Local, immutable, Link Length
    theta = zeros(n, 1); % Rotation angle at each node
    theta_p = zeros(n, 1); % First derivative of theta
    theta_p_old = ones(n, 1) * 1e-6; % Previous iteration of theta_p
    r = zeros(n, 2); % Position vector at each node

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
    m = M * magMoments; % Scale the provided magnetic moments by M

    % Create the magnetic field template
    xLow = -0.1;   yLow = -0.1;
    xHigh = 0.1;   yHigh = 0.1;
    xIncr = 0.005;  yIncr = 0.005; % Step size

    % Create the meshgrid used for magnetic simulation
    [Fieldx, Fieldy] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh);

    % Define the magnetic field Bx and By over the grid
    % Use Bxin and Byin as constants or define functions for non-uniform fields
    Bx = Bxin * ones(size(Fieldx)); % Replace with desired field distribution
    By = Byin * ones(size(Fieldy)); % Replace with desired field distribution

    % Loop parameters
    loops = 1000; % Maximum number of iterations
    mu = 0.5; % Numerical damping parameter

    R = zeros(2, 2, n); % Rotation matrices at each node
    R_hat = zeros(2, 2, n); % Derivative of rotation matrices wrt theta

    r(1,:) = [0, 0]; % Starting position at the origin

    for j = 1:loops
        % Compute rotation matrices and their derivatives
        for i = 1:n
            % Rotation Matrix for current theta
            cos_theta = cos(theta(i));
            sin_theta = sin(theta(i));
            R(:, :, i) = [cos_theta, -sin_theta;
                          sin_theta,  cos_theta];
            m_def(i, :) = (R(:, :, i) * m(i, :)')';

            % Derivative of Rotation Matrix wrt theta
            R_hat(:, :, i) = [-sin_theta, -cos_theta;
                               cos_theta, -sin_theta];
            m_hat(i, :) = (R_hat(:, :, i) * m(i, :)')';
        end

        % Map to Cartesian coordinates
        dx = delta_s(1) * cos(theta);
        dy = delta_s(1) * sin(theta);
        r(1, :) = [0, 0];
        r(2:end, 1) = cumsum(dx(1:end-1)) + r(1,1);
        r(2:end, 2) = cumsum(dy(1:end-1)) + r(1,2);

        % Interpolate Bx and By at positions r(:,1), r(:,2)
        % Handle points outside the interpolation grid by setting them to zero
        Bx_interp = interp2(Fieldx, Fieldy, Bx, real(r(:,1)), real(r(:,2)), 'linear', 0);
        By_interp = interp2(Fieldx, Fieldy, By, real(r(:,1)), real(r(:,2)), 'linear', 0);
        B_interp = [Bx_interp, By_interp];

        % Compute the expression for updating theta_p
        EI = E * I; % Bending stiffness
        gravity_term = (rho * g * A / EI) * cos(theta);
        magnetic_term = (theta_p_old.^2) .* (sum(m_hat .* B_interp, 2) / EI);
        expr = gravity_term + magnetic_term;

        % Update theta_p with numerical damping
        theta_p = (1 - mu) * theta_p_old + mu * cbrt(expr);

        % Integrate theta_p to get theta
        theta(2:end) = theta(1:end-1) + ((theta_p(2:end) + theta_p(1:end-1)) / 2) * delta_s(1);

        % Convergence Clause
        err_tol = 1e-6;
        if norm(theta_p - theta_p_old) < err_tol
            fprintf('Converged at iteration %d with error norm %e\n', j, norm(theta_p - theta_p_old));
            break;
        end

        if j == loops
            fprintf('Maximum iterations reached without convergence. Error norm: %e\n', norm(theta_p - theta_p_old));
        end

        % Update theta_p_old for the next iteration
        theta_p_old = theta_p;
    end

    % Return the x and y positions
    xpos = real(r(:, 1));
    ypos = real(r(:, 2));
end

% --- Helper Functions ---
function y = cbrt(x)
    y = sign(x) .* abs(x).^(1/3);
end
