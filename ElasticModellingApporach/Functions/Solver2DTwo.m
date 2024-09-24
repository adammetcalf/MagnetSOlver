function [xpos, zpos] = Solver2DTwo(n, L, magFieldFunc, magMoments)
    % Numeric
    s = linspace(0, L, n)'; % Column vector
    delta_s = repmat(s(2) - s(1), n, 1); % Local, immutable, link length
    theta = zeros(n, 1); % Rotation angle at each node
    theta_p = zeros(n, 1); % First derivative of theta
    theta_p_old = ones(n, 1) * 1e-6; % Previous iteration of theta_p
    r = zeros(n, 2); % Position vector at each node

    % Geometric
    radius = 2e-3; % Radius [m]
    A = pi * radius^2; % Cross-sectional area [m^2]
    V = A; % Volume per unit length [m^2]
    I = 0.25 * pi * radius^4; % Second moment of area [m^4]

    % Material
    E = 10e5; % Young's modulus [N/m^2]
    rho = 1100; % Density [kg/m^3]
    g = -9.81; % Acceleration due to gravity [m/s^2]

    EI = E*I;

    % Magnetic
    mu0 = 4 * pi * 1e-7; % Vacuum permeability (H/m or N/A^2)
    remanence = 0.10; % Saturation remanent magnetic field [T]
    M = remanence * V / mu0; % Magnitude of magnetic moment (A·m^2)
    m = zeros(n, 2); % Spatial vector of referential magnetic moment
    m_def = zeros(n, 2); % Spatial vector of deformed magnetic moment
    m_hat = zeros(n, 2); % Derivative wrt theta of deformed magnetic moment

    % Initialize magnetic moments
    m = M * magMoments; % Scale the provided magnetic moments by M

    % Loop parameters
    loops = 1000; % Maximum number of iterations
    mu = 0.5; % Numerical damping parameter

    R = zeros(2, 2, n); % Rotation matrices at each node
    R_hat = zeros(2, 2, n); % Derivative of rotation matrices wrt theta

    r(1,:) = [0, 0]; % Starting position at the origin

    for j = 1:loops
        % Compute rotation matrices and their derivatives
        for i = 1:n
            % Rotation matrix for current theta
            cos_theta = cos(theta(i));
            sin_theta = sin(theta(i));
            R(:, :, i) = [cos_theta, -sin_theta;
                          sin_theta,  cos_theta];
            m_def(i, :) = (R(:, :, i) * m(i, :)')';

            % Derivative of rotation matrix wrt theta
            R_hat(:, :, i) = [-sin_theta, -cos_theta;
                               cos_theta, -sin_theta];
            m_hat(i, :) = (R_hat(:, :, i) * m(i, :)')';
        end

         % Map to Cartesian coordinates
        dx = delta_s(1) * cos(theta);
        dz = delta_s(1) * sin(theta);
        r(1, :) = [0, 0];
        r(2:end, 1) = cumsum(dx(1:end-1)) + r(1,1);
        r(2:end, 2) = cumsum(dz(1:end-1)) + r(1,2);
    
        % Compute magnetic field at current positions
        [Bx_interp, Bz_interp] = magFieldFunc(r(:,1), r(:,2));
        B_interp = [Bx_interp(:), Bz_interp(:)];

        % Compute the expression for updating theta_p
        gravity_term = (rho * g * A / EI) * cos(theta);
        magnetic_term = (sum(m_hat .* B_interp, 2) / EI);
        expr = gravity_term + magnetic_term;

        % Update theta_p with numerical damping
        theta_p = (1 - mu) * theta_p_old + mu * cbrt(expr);

        % Integrate theta_p to get theta
        theta(2:end) = theta(1:end-1) + ((theta_p(2:end) + theta_p(1:end-1)) / 2) .* delta_s(1);

        % Convergence check
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
    zpos = real(r(:, 2));
end

% --- Helper Functions ---
function y = cbrt(x)
    y = sign(x) .* abs(x).^(1/3);
end