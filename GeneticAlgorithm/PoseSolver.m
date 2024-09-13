classdef PoseSolver
    %POSESOLVER Optimizer to get the pose under current conditions
    %   This solver uses optimization to solve for the pose of the tentacle
    %   under current magnetic conditions

    properties (Access = public)
        World World; % A world object
        lowerBounds double;
        upperBounds double;
        initialAngles double;
        optimizeAngleIndices;
    end

    methods

        function obj = PoseSolver(World, optimizeAngleIndices)
            %POSESOLVER Construct an instance of this class

            % Inject the world
            obj.World = World;

            initialAngles = obj.World.getJointAngles();
            obj.initialAngles = initialAngles;

            if nargin < 2
                % By default, optimize over all angles
                obj.optimizeAngleIndices = 1:size(initialAngles, 2);
            else
                obj.optimizeAngleIndices = optimizeAngleIndices;
            end

            % Get the number of joints and angles
            numJoints = size(initialAngles, 1);
            numAnglesToOptimize = length(obj.optimizeAngleIndices);

            % Set the lower and upper bounds for the angles to optimize
            obj.lowerBounds = repmat(-180, numJoints, numAnglesToOptimize);
            obj.upperBounds = repmat(180, numJoints, numAnglesToOptimize);
        end

        function newAngles = optimizeJoints(obj)
            initialAngles = obj.World.getJointAngles(); % size numJoints x numAnglesPerJoint

            % Extract the angles to optimize
            anglesToOptimize = initialAngles(:, obj.optimizeAngleIndices);

            % Flatten the angles into a vector
            initialGuess = anglesToOptimize(:);

            options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'none');

            % Flatten lower and upper bounds
            lb = obj.lowerBounds(:);
            ub = obj.upperBounds(:);

            [optimizedAngles, ~] = fmincon(@obj.objectiveFunction, initialGuess, [], [], [], [], lb, ub, [], options);

            % Reshape optimizedAngles back into numJoints x numAnglesToOptimize
            optimizedAngles = reshape(optimizedAngles, size(anglesToOptimize));

            % Form the newAngles matrix
            newAngles = obj.initialAngles; % Start with initial angles
            newAngles(:, obj.optimizeAngleIndices) = optimizedAngles; % Replace the optimized angles
        end

        function cost = objectiveFunction(obj, angles)
            % angles is a vector of the angles being optimized, reshape accordingly
            numJoints = size(obj.initialAngles, 1);
            numAnglesToOptimize = length(obj.optimizeAngleIndices);
            angles = reshape(angles, numJoints, numAnglesToOptimize);

            % Form the complete angles matrix
            completeAngles = obj.initialAngles; % Start with initial angles
            completeAngles(:, obj.optimizeAngleIndices) = angles; % Replace the optimized angles

            % Update the angles in the World
            obj.World = obj.World.UpDateAngles(completeAngles);  % Update the angles

            forcesTorques = obj.World.getForcesTorques();
            forcesTorques = abs(forcesTorques);
            cost = sum(((forcesTorques .* 100)), "all");
        end

    end
end