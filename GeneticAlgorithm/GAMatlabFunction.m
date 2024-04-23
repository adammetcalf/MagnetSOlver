classdef GAMatlabFunction
    %GAMATLABFUNCTION Summary of this class goes here
    %   This solver uses a genetic algorithm to solve for the pose of the
    %   tentacle under current magnetic conditions
    
    properties 
        World                               % A world object
        lowerBounds double = -180;
        upperBounds double = 180;
        initialThetas double;
    end
    
    methods
        function obj = GAMatlabFunction(World)
            %GAMATLABFUNCTION Construct an instance of this class
            
            % Inject the world
            obj.World = World;
 
            initialAngles = obj.World.getJointAngles();
            obj.initialThetas = initialAngles(:,1); % Currently optimising only the alpha angles
            % #TODO optimise the theta angles too

            %Get the number of rows in angles
            numAngles = size(initialAngles, 1);

            % Set the lower and upper bounds for the theta angles
            obj.lowerBounds = repmat(obj.lowerBounds, numAngles, 1);
            obj.upperBounds = repmat(obj.upperBounds, numAngles, 1);
        end
        
        function newAngles = optimizeJoints(obj)
            initialAngles = obj.World.getJointAngles();
            initialAngles = initialAngles(:,2); % Starting angles for optimization
            nVars = length(initialAngles);  % Number of variables to optimize
            
            % Genetic Algorithm options setup
            options = optimoptions('ga', ...
                                   'PopulationSize', 100, ...
                                   'MaxGenerations', 100, ...
                                   'CrossoverFcn', @crossoverintermediate, ...
                                   'MutationFcn', @mutationadaptfeasible, ...
                                   'Display', 'iter');  % Iterative display of the process

            % Running the Genetic Algorithm
            [optimizedAngles, ~] = ga(@obj.objectiveFunctionGA, nVars, [], [], [], [], obj.lowerBounds, obj.upperBounds, [], options);
            newAngles = [obj.initialThetas, optimizedAngles'];  % Combine the constant and optimized angles
        end
        
        function cost = objectiveFunctionGA(obj, angles)
            % This function forms the complete angle set and computes the cost
            angles = [obj.initialThetas, angles']; % Form complete angles
            obj.World = obj.World.UpDateAngles(angles);  % Update the angles
            forcesTorques = obj.World.getForcesTorques(); 
            forcesTorques = abs(forcesTorques);
            cost = sum((forcesTorques.*100), "all");
        end
    end
end

