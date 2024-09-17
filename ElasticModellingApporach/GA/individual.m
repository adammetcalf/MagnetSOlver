classdef individual
    %individual class. An indivudual is a single member of a population

    properties (Access = private)
        Field double;
        Fitness double;
    end

    methods (Access = public)

        %% Constructor
        function obj = individual(world)

            % constructs an instance of individual, and injects the world
            % object
            obj.Field = world.getField();               % init the Field variable
            obj = randomise(obj);                       % randomise the filed (initial guess)
            obj = determineFitness(obj, world);         % measure the fitness of this solution

        end

        %% Accessors
        function Field = getField(obj)
            Field = obj.Field;
        end

        function Fitness = getFitness(obj)
            Fitness = obj.Fitness;
        end

        function obj = updateField(obj,Field,world)

            obj.Field = Field;                          % Update Field
            obj = determineFitness(obj, world);         % measure the fitness of this solution
        end


    end

    methods (Access = private)

        function obj = randomise(obj)

            % Define limits
            lower_limit = -0.0025;  
            upper_limit = 0.0025; 
            
            % Randomize a and b between the limits
            a = lower_limit + (upper_limit - lower_limit) * rand;
            b = lower_limit + (upper_limit - lower_limit) * rand;
            
            % Create the vector x
            obj.Field = [a, b];
        end


        function obj = determineFitness(obj, world)

            % inject new randomised field into world
            world = world.updateField(obj.Field);

            % obtain resultant positions
            [DesiredXPos,DesiredYPos,xPos,yPos] = world.getPositions();

            % determine fintess
            obj.Fitness = ((xPos - DesiredXPos).^2 + (yPos - DesiredYPos).^2);

        end

    end

end