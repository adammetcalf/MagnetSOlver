classdef PoseSolver
    %POSESOLVER Optimiser to get the pose under current conditions
    %   This solver uses optimisation to solve for the pose of the tentacle
    %   under current magnetic conditions
    
    properties (Access = public)
        World World; % A world object
        lowerBounds double;
        upperBounds double;
        initialThetas double;
    end
    
    methods

        function obj = PoseSolver(World)
            %POSESOLVER Construct an instance of this class
            
            obj.World = World;
 
            initialAngles = obj.World.getJointAngles();
            obj.initialThetas = initialAngles(:,1); % Currently optimising only the alpha angles

        end
        
        function newAngles = optimizeJoints(obj)
            initialAngles = obj.World.getJointAngles();
            initialAngles = initialAngles(:,2);
            options = optimoptions('fmincon', 'Algorithm', 'sqp');

            %[optimizedAngles,~] = fmincon(@obj.objectiveFunction, initialAngles, [], [], [], [], obj.lowerBounds, obj.upperBounds, @obj.constraints, options);
            [optimizedAngles,~] = fmincon(@obj.objectiveFunction, initialAngles, [], [], [], [], [], [], [], options);
            newAngles = [obj.initialThetas,optimizedAngles];
        end
        
        function cost = objectiveFunction(obj, angles)
            angles=[obj.initialThetas,angles]; %form complete angles
            obj.World = obj.World.UpDateAngles(angles);  % Update the angles
            forcesTorques = obj.World.getForcesTorques(); 
            forcesTorques = abs(forcesTorques);
            cost = sum(forcesTorques,"all"); % Example cost function
        end

    end
end
