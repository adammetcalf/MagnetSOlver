classdef TentacleBuilder
    %TENTACLEBUILDER Builds a tentacle for the simulation
    
    properties
        tentacle Tentacle;  % Holds a tentacle
    end
    
    methods
        function obj = TentacleBuilder(inputArg1,inputArg2)
            %TENTACLEBUILDER Launch UI
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

