classdef PoseSolver
    %POSESOLVER Optimiser to get the pose under current conditions
    %   This solver uses optimisation to solve for the pose of the tentacle
    %   under current magnetic conditions
    
    properties (Access = public)
        World World; % A world object
        lowerBounds double = -180;
        upperBounds double = 180;
        initialThetas double;
    end
    
    methods

        function obj = PoseSolver(World)
            %POSESOLVER Construct an instance of this class
            
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
            initialAngles = initialAngles(:,2);
            options = optimoptions('fmincon', 'Algorithm', 'sqp','Display', 'none'); % Display: none suppresses all the annoying messages

            [optimizedAngles,~] = fmincon(@obj.objectiveFunction, initialAngles, [], [], [], [], obj.lowerBounds, obj.upperBounds, [], options);
            newAngles = [obj.initialThetas,optimizedAngles];
        end
        
        function cost = objectiveFunction(obj, angles)
            angles=[obj.initialThetas,angles]; %form complete angles
            obj.World = obj.World.UpDateAngles(angles);  % Update the angles
            forcesTorques = obj.World.getForcesTorques(); 
            forcesTorques = abs(forcesTorques);
            cost = sum(((forcesTorques.*100)),"all");
            
        end

    end
end
