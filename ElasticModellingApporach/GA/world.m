classdef world
% The world contains the magnetically actuated elastica and the Actuating
% magnetic field

    %% private properties
    properties (Access = private)
        DesiredXPos double;        % A vector of the desired x-positions
        DesiredYPos double;        % A vector of the desired y-positions
        xPos double;        % A vector of the x-positions
        yPos double;        % A vector of the y positions
        magMoments double;  % A matrix of the local magnetic moment directions
        Field double;       % A vector containg the x, y field (this is the component that changes to solve the GA)
        L;                  % total elastica length

        %Magnetic field template
        Fieldx double;      % Holds part of the magnetic filed
        Fieldy double;      % Holds part of the magnetic field

        % MAgnetic field
        Bx double;          % The magneitc field
        By double;          % The magnetic field

        %Geometric properties
        radius = 2e-3; % Radius [m]
        A double; % Cross Sectional Area [m^2]
        V double; % Volume per unit length [m^2]
        I double;  % Second moment of area [m^4]

        % Material
        E = 10e5; % Young's modulus [N/m^2]
        rho = 1100; % Density [kg/m^3]
        g = -9.81; % Acceleration due to gravity [m/s^2]

        % Magnetic
        mu0 = 4 * pi * 1e-7; % Air magnetic permeability (H/m or N/A^2)
        remanence = 0.10; % Saturation Remanent Magnetic Field [T]
        M double; % Magnitude of magnetic moment (A·m^2)
    end

    %% Public methods
    methods (Access = public)

        %% Constructor
        function obj = world(DesiredXPos,DesiredYPos,magMoments, Length)

            % init variables
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

            % Create the magnetic field
            obj = createField(obj);
        end

        %% Accessors
        function Field = getField(obj)
            Field = obj.Field;
        end

        function [DesiredXPos,DesiredYPos,xPos,yPos] = getPositions(obj)
            xPos = obj.xPos;
            yPos = obj.yPos;
            DesiredXPos = obj.DesiredXPos;
            DesiredYPos = obj.DesiredYPos;
        end

        function obj = updateField(obj,Field)
        % Update the field and evaluate the world.
            obj.Field = Field;
            obj = evaluateWorld(obj);
        end

    end

    %% private methods
    methods (Access = private)
      
        function obj = createField(obj)

            % Create the magnetic field template
            xLow = -0.1;   yLow = -0.1;
            xHigh = 0.1;   yHigh = 0.1;
            xIncr = 0.005;  yIncr = 0.005; % Step size

            % Create the meshgrid used for magnetic simulation
            [obj.Fieldx, obj.Fieldy] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh);
            
            % Define the magnetic field Bx and By over the grid
            obj.Bx = obj.Field(1) * ones(size(obj.Fieldx)); 
            obj.By = obj.Field(2) * ones(size(obj.Fieldy)); 

        end

        % Evaluates the response of the tentacle to the Field
        function obj = evaluateWorld(obj)

            % Number of nodes in elastica
            n = length(obj.xPos);

            s = linspace(0, obj.L, n)'; % Column vector
            delta_s = repmat(s(2) - s(1), n, 1); % Local, immutable, Link Length
            theta = zeros(n, 1); % Rotation Matrix from d3 to e3
            theta_p = zeros(n, 1); % Differential of R wrt s
            theta_p_old = ones(n, 1) * 1e-6; % Previous iteration of theta_p
            r = zeros(n, 2); % Vector position as a function of "s"

            % Initialize magnetic moments
            m_def = zeros(n, 2); % Spatial vector of deformed magnetic moment
            m_hat = zeros(n, 2); % Differential wrt theta of deformed magnetic moments
            m = obj.M * obj.magMoments; % Scale the provided magnetic moments by M

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
                Bx_interp = interp2(obj.Fieldx, obj.Fieldy, obj.Bx, real(r(:,1)), real(r(:,2)), 'linear', 0);
                By_interp = interp2(obj.Fieldx, obj.Fieldy, obj.By, real(r(:,1)), real(r(:,2)), 'linear', 0);
                B_interp = [Bx_interp, By_interp];
        
                % Compute the expression for updating theta_p
                EI = obj.E * obj.I; % Bending stiffness
                gravity_term = (obj.rho * obj.g * obj.A / EI) * cos(theta);
                magnetic_term = (theta_p_old.^2) .* (sum(m_hat .* B_interp, 2) / EI);
                expr = gravity_term + magnetic_term;
        
                % Update theta_p with numerical damping
                theta_p = (1 - mu) * theta_p_old + mu * cbrt(expr);
        
                % Integrate theta_p to get theta
                theta(2:end) = theta(1:end-1) + ((theta_p(2:end) + theta_p(1:end-1)) / 2) * delta_s(1);

                % Convergence Clause
                err_tol = 1e-6;
                if norm(theta_p - theta_p_old) < err_tol
                    %fprintf('Converged at iteration %d with error norm %e\n', j, norm(theta_p - theta_p_old));
                    break;
                end
        
                if j == loops
                    %fprintf('Maximum iterations reached without convergence. Error norm: %e\n', norm(theta_p - theta_p_old));
                end
        
                % Update theta_p_old for the next iteration
                theta_p_old = theta_p;
            end

            % Return the x and y positions
            obj.xPos = real(r(:, 1));
            obj.yPos = real(r(:, 2));

        end

        % Helper function (prevents complex numbers)
        function y = cbrt(x)
            y = sign(x) .* abs(x).^(1/3);
        end
    end




end