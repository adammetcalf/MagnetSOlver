classdef PoseSolver
    %POSESOLVER Optimiser to get the pose under current conditions
    %   This solver uses optimisation to solve for the pose of the tentacle
    %   under current magnetic conditions
    
    properties (Access = public)
        World World; % A world object
        lowerBounds double = -180;
        upperBounds double = 180;
    end
    
    methods

        function obj = PoseSolver(World)
            %POSESOLVER Construct an instance of this class
            
            obj.World = World;
 
            initialAngles = obj.World.getJointAngles();
            % Assuming 'initialAngles' is a Nx2 matrix, where N is the number of joints
            
            obj.lowerBounds = -180 * ones(size(initialAngles)); % Lower bounds of -180 degrees for each angle
            obj.upperBounds = 180 * ones(size(initialAngles)); % Upper bounds of 180 degrees for each angle

        end
        
        function optimizedAngles = optimizeJoints(obj)
            initialAngles = obj.World.getJointAngles();
            options = optimoptions('fmincon', 'Algorithm', 'sqp');
            [optimizedAngles,~] = fmincon(@obj.objectiveFunction, initialAngles, [], [], [], [], obj.lowerBounds, obj.upperBounds, @obj.constraints, options);
            obj = obj.World.UpDateAngles(optimizedAngles);
        end
        
        function cost = objectiveFunction(obj, angles)
            obj = obj.World.UpDateAngles(angles);  % Update the angles
            forcesTorques = obj.World.getForcesTorques(); % You need to implement this
            cost = sum(abs(forcesTorques),all) % Example cost function
        end
        
        function [c, ceq] = constraints(obj, angles)
            % No equality constraints
            ceq = [];
            
            % Inequality constraints
            upperLimit = 30; % degrees, convert to radians if using radians in your system
            lowerLimit = -30; % degrees, convert to radians if using radians in your system
            
            % Constraints for the second column to be within +/- 30 degrees
            c_upper = angles(:, 2) - upperLimit; % Should be less than or equal to 0
            c_lower = lowerLimit - angles(:, 2); % Should be less than or equal to 0
            
            % Combine constraints
            c = [c_upper; c_lower];
        end

    end
end

%function optimized_sum = optimizeMatrixA(A_initial, B_function)
    % A_initial is the initial state of matrix A [N x 2]
    % B_function is a function handle that computes matrix B and its sum based on A

    % Objective function: computes the sum of elements in B based on A
   % objectiveFunction = @(A) computeSumB(A, B_function);
    
    % Initial guess for the second column of A
   % A0 = A_initial(:, 2);

    % Linear inequality constraints (empty for this problem)
   % A_eq = [];
   % b_eq = [];
    
    % Bounds for the second column of A, ensuring values are within +/- 30 of the initial values
  %  lb = A_initial(:, 2) - 30;
  %  ub = A_initial(:, 2) + 30;

    % Linear equality constraints (ensure the first column of A is unchanged)
    % We are not modifying the first column, so we can ignore this in the optimization

    % Options for the optimization (optional)
  %  options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    
    % Solve the optimization problem
  %  [A_optimized, fval] = fmincon(objectiveFunction, A0, [], [], A_eq, b_eq, lb, ub, [], options);
    
    % Update A with optimized values
  %  A_final = A_initial;
  %  A_final(:, 2) = A_optimized;

    % Return the optimized sum of B
   % optimized_sum = fval;
%end

%function sumB = computeSumB(A, B_function)
%    % Compute matrix B based on A
%    B = B_function(A);
    % Compute the sum of elements in B
%    sumB = sum(B, 'all');
%end
