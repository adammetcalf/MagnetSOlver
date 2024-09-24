classdef world
% The world contains the magnetically actuated elastica and the actuating
% magnetic field.

    %% Private Properties
    properties (Access = private)
        DesiredXPos double;   % A vector of the desired x-positions
        DesiredYPos double;   % A vector of the desired y-positions
        xPos double;          % A vector of the x-positions
        yPos double;          % A vector of the y positions
        magMoments double;    % A matrix of the local magnetic moment directions
        Field double;         % A vector containing the x, y field
        L;                    % Total elastica length
        orientation;          % Starting orientation (0 for horizontal, 1 for vertically downwards)

        % Magnetic field template
        Fieldx double;        % X-coordinates for magnetic field grid
        Fieldy double;        % Y-coordinates for magnetic field grid

        % Magnetic field components
        Bx double;            % X-component of the magnetic field
        By double;            % Y-component of the magnetic field

        % Geometric properties
        radius = 2e-3;        % Radius [m]
        A double;             % Cross-sectional area [m^2]
        V double;             % Volume per unit length [m^2]
        I double;             % Second moment of area [m^4]

        % Material properties
        E = 10e5;             % Young's modulus [N/m^2]
        rho = 1100;           % Density [kg/m^3]
        g = -9.81;            % Acceleration due to gravity [m/s^2]

        % Magnetic properties
        mu0 = 4 * pi * 1e-7;  % Air magnetic permeability [H/m]
        remanence = 0.10;     % Saturation remanent magnetic field [T]
        M double;             % Magnitude of magnetic moment [A·m^2]
    end

    %% Public Methods
    methods (Access = public)

        %% Constructor
        function obj = world(DesiredXPos, DesiredYPos, magMoments, Length, orientation)
            % Initialize variables
            obj.xPos = DesiredXPos;
            obj.yPos = DesiredYPos;
            obj.DesiredXPos = DesiredXPos;
            obj.DesiredYPos = DesiredYPos;
            obj.magMoments = magMoments;
            obj.L = Length;
            obj.A = pi * obj.radius^2;
            obj.V = obj.A;
            obj.I = 0.25 * pi * obj.radius^4;
            obj.Field = [0.0, 0.0];
            obj.M = obj.remanence * obj.V / obj.mu0;
            obj.orientation = orientation;

            % Create the magnetic field
            obj = createField(obj);
        end

        %% Accessors
        function Field = getField(obj)
            Field = obj.Field;
        end

        function [DesiredXPos, DesiredYPos, xPos, yPos] = getPositions(obj)
            xPos = obj.xPos;
            yPos = obj.yPos;
            DesiredXPos = obj.DesiredXPos;
            DesiredYPos = obj.DesiredYPos;
        end

        function [magMoments, Length] = getMagLength(obj)
            magMoments = obj.magMoments;
            Length = obj.L;
        end

        function obj = updateField(obj, Field)
            % Update the field and evaluate the world
            obj.Field = Field;

            % Update Bx and By based on the new Field
            obj.Bx = obj.Field(1) * ones(size(obj.Fieldx));
            obj.By = obj.Field(2) * ones(size(obj.Fieldy));

            obj = evaluateWorld(obj);
        end

        function orientation = getOri(obj)

            orientation = obj.orientation;
        end

    end

    %% Private Methods
    methods (Access = private)

        function obj = createField(obj)
            % Create the magnetic field template
            xLow = -0.1;   yLow = -0.1;
            xHigh = 0.1;   yHigh = 0.1;
            xIncr = 0.005; yIncr = 0.005; % Step size

            % Create the meshgrid used for magnetic simulation
            [obj.Fieldx, obj.Fieldy] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh);

            % Define the magnetic field Bx and By over the grid
            obj.Bx = obj.Field(1) * ones(size(obj.Fieldx));
            obj.By = obj.Field(2) * ones(size(obj.Fieldy));
        end

        % Evaluates the response of the tentacle to the field
        function obj = evaluateWorld(obj)

            % Number of nodes in elastica
            n = length(obj.xPos);

            s = linspace(0, obj.L, n)'; % Column vector
            delta_s = repmat(s(2) - s(1), n, 1); % Local, immutable link length
            theta = zeros(n, 1); % Rotation angle at each node
            theta_p = zeros(n, 1); % First derivative of theta
            theta_p_old = ones(n, 1) * 1e-6; % Previous iteration of theta_p
            r = zeros(n, 2); % Position vector at each node

            % Initialize magnetic moments
            m_def = zeros(n, 2); % Spatial vector of deformed magnetic moment
            m_hat = zeros(n, 2); % Differential wrt theta of deformed magnetic moments
            m = obj.M * obj.magMoments; % Scale the provided magnetic moments by M

            % Loop parameters
            loops = 100; % Maximum number of iterations
            mu = 0.5;     % Numerical damping parameter

            R = zeros(2, 2, n);     % Rotation matrices at each node
            R_hat = zeros(2, 2, n); % Derivative of rotation matrices wrt theta

            r(1, :) = [0, 0]; % Starting position at the origin

            % Set initial orientation based on the obj.orientation property
            if obj.orientation == 0
                % Horizontal orientation
                theta = zeros(n, 1); % Tentacle starts horizontally along the x-axis
            elseif obj.orientation == 1
                % Vertically downwards orientation
                theta = (-pi / 2) * ones(n, 1); % Tentacle hangs downwards along the negative y-axis
            else
                error('Invalid orientation. Use 0 for horizontal, 1 for vertically downwards.');
            end

            for j = 1:loops

                % Compute rotation matrices and their derivatives
                cos_theta = cos(theta);
                sin_theta = sin(theta);

                % Initialize R and R_hat as 2 x 2 x n arrays
                R = zeros(2, 2, n);
                R_hat = zeros(2, 2, n);

                % Populate R
                R(1, 1, :) = cos_theta;
                R(1, 2, :) = -sin_theta;
                R(2, 1, :) = sin_theta;
                R(2, 2, :) = cos_theta;

                % Populate R_hat
                R_hat(1, 1, :) = -sin_theta;
                R_hat(1, 2, :) = -cos_theta;
                R_hat(2, 1, :) = cos_theta;
                R_hat(2, 2, :) = -sin_theta;

                % Compute m_def and m_hat
                % m is n x 2, we need to perform matrix multiplication for each slice
                for i = 1:n
                    m_def(i, :) = (R(:, :, i) * m(i, :)')';
                    m_hat(i, :) = (R_hat(:, :, i) * m(i, :)')';
                end

                % Map to Cartesian coordinates
                dx = delta_s(1) * cos_theta;
                dy = delta_s(1) * sin_theta;
                r(1, :) = [0, 0];
                r(2:end, 1) = cumsum(dx(1:end-1)) + r(1, 1);
                r(2:end, 2) = cumsum(dy(1:end-1)) + r(1, 2);

                % Interpolate Bx and By at positions r(:,1), r(:,2)
                Bx_interp = interp2(obj.Fieldx, obj.Fieldy, obj.Bx, real(r(:, 1)), real(r(:, 2)), 'linear', 0);
                By_interp = interp2(obj.Fieldx, obj.Fieldy, obj.By, real(r(:, 1)), real(r(:, 2)), 'linear', 0);
                B_interp = [Bx_interp, By_interp];

                % Compute the expression for updating theta_p
                EI = obj.E * obj.I; % Bending stiffness
                gravity_term = (obj.rho * obj.g * obj.A / EI) .* cos_theta;
                magnetic_term = (theta_p_old.^2) .* (sum(m_hat .* B_interp, 2) / EI);
                expr = gravity_term + magnetic_term;

                % Update theta_p with numerical damping
                theta_p = (1 - mu) * theta_p_old + mu * obj.cbrt(expr);

                % Integrate theta_p to get theta
                theta(2:end) = theta(1:end-1) + ((theta_p(2:end) + theta_p(1:end-1)) / 2) * delta_s(1);

                % Convergence Clause
                err_tol = 1e-6;
                if norm(theta_p - theta_p_old) < err_tol
                    break;
                end

                % Update theta_p_old for the next iteration
                theta_p_old = theta_p;
            end

            % Update the x and y positions
            obj.xPos = real(r(:, 1));
            obj.yPos = real(r(:, 2));

        end

        % Helper function to compute the cube root, handling negative values
        function y = cbrt(obj, x)
            y = sign(x) .* abs(x).^(1/3);
        end

    end

end
