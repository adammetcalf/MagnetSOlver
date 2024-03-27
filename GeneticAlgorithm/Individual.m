classdef Individual
    %individual class. An individual is a single member of a population
    
    properties (Access = private)
        Angles double;
        fitness double;
    end
    
    methods (Access = public)

        %% Constructor
        function obj = Individual(World)
            
            %INDIVIDUAL Construct an instance of this class and inject a
            %world object
           
            obj.Angles = World.getJointAngles();  %init Angles variable
            obj = randomise(obj); %randomise the joint angles (initial guess)
            obj = determineFitness(obj, World); % measure fitness of this solution
        end

        function obj = updateAngles(obj, angles, World)

            %Update joint angles
            obj.Angles = angles;
            
            obj = determineFitness(obj, World); % measure fitness of this solution
        end

        %% Accessors
        function Angles = getAngles(obj)
            Angles = obj.Angles;
        end

        function fitness = getFitness(obj)
            fitness = obj.fitness;
        end

    end
       
    methods (Access = private)

        function obj = randomise(obj)         
  
            % Randomise all the alpha values to between -180 and 180.
            obj.Angles(:,2) = -180 + (360)*rand(size(obj.Angles,1), 1);
        end

        function obj = determineFitness(obj, World)

            % Update joint angles in world
            World = World.UpDateAngles(obj.Angles);

            % Get Forces and Torques
            ForcesTorques = World.getForcesTorques();

            % Evaluate Fitness
            ForcesTorques = abs(ForcesTorques);
            obj.fitness = sum(ForcesTorques,"all");
        end
    end
end

