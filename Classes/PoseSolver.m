classdef PoseSolver
    %POSESOLVER Optimiser to get the pose under current conditions
    %   This solver uses optimisation to solve for the pose of the tentacle
    %   under current magnetic conditions
    
    properties (Access = public)
        World World; % A world object
    end
    
    methods
        function obj = PoseSolver(inputArg1,inputArg2)
            %POSESOLVER Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

