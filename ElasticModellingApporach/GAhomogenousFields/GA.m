classdef GA
% This runs a Genetic algorithm to solve for the field fo r a desired
% tentacle shape.

    %% private properties
    properties (Access = private)
        world world;
        bestIndividual individual;
        NextBestIndividual individual;

    end

    %% public methods
    methods (Access = public)

        %% Constructor
        function obj = GA(world)
            
            obj.world = world;
            obj.bestIndividual = individual(world);
            obj.NextBestIndividual = individual(world);


        end

        function BestField = Train(obj,opt)

            %init epoch count
            epoch = 0;        
            
            % init angle change count
            FieldNoChange = 0;

            if opt %if optimising, use less epochs (2)
                endEpoch = 1;
            else
                endEpoch =5;
            end

            while ((epoch < endEpoch) && (FieldNoChange <= 1))
                disp("Epoch: "+ num2str(epoch));

                % Best angles at start of epoch
                FieldStart = obj.bestIndividual.getField();

                % perform 1 epoch of training
                obj = obj.performEpoch();

                if opt
                % Optimise best individual
                    obj = obj.Optimise();
                end

                % increment epoch count
                epoch = epoch+1;

                % Best angles at end of epoch
                FieldEnd = obj.bestIndividual.getField();

                % Compare best angles at start and end of epoch
                if FieldStart == FieldEnd
                    FieldNoChange = FieldNoChange+1;
                else
                    FieldNoChange =0;
                end

                disp(FieldNoChange);

            end
            
            % Ensure that best angles have propogated through (probably
            % unnecessary)
            BestField = obj.bestIndividual.getField();
        end


    end


    %% Private methods
    methods (Access = private)

        % perform the epoch
        function obj = performEpoch(obj)
            
            % initialise epoch
            pop = population(obj.world);
            
            %inject best individual from previous epoch
            pop = pop.injectBest(obj.bestIndividual, obj.NextBestIndividual);

            % Perform the evolution cycles
            obj = obj.performEvolution(pop);
 
        end

        % perform the evolution until there is no further convolution
        function obj = performEvolution(obj, pop)

            % init convolution check count, which is used to track
            % improvements across evloution cycles.

            Evolution = 0;

            convCheckCount = 0;
            while convCheckCount < 5  
                

                % Get the best individual fitnesss
                [~, fitness, ~] = pop.getBest();

                % Evolve the population
                pop = pop.Evolve(obj.world, Evolution);

                % Get the best individual after evolution
                [~, fitness2, ~] = pop.getBest();

                % if there has been an improvement after evolution, reset
                % convolution count to zero, else increment
                if abs(fitness2 - fitness) < 1e-8
                    convCheckCount = convCheckCount + 1;
                else
                    convCheckCount = 0;
                end
                Evolution = Evolution+1;


                disp("Evolution: "+ num2str(Evolution)+ ", Fitness: "+ num2str(fitness2));
            end 

            % End of the evolution cycles for this epoch (ie no improvement
            % for 20 consequtive evolutions). 

            %update the best individual
            popBest = pop.getBestIndividual();
            popFit = popBest.getFitness();
            oldFit = obj.bestIndividual.getFitness();

            if popFit < oldFit
                obj.bestIndividual = pop.getBestIndividual();
            end

        end

        function obj = Optimise(obj)
            % Extract the current best field
            initial_field = obj.bestIndividual.getField();
        
            % Desired positions from the user interface or desired shape
            [x_positionsUI, y_positionsUI, ~, ~] = obj.world.getPositions();
        
            % Get magnet moments and length
            [magn_momentsUI, L] = obj.world.getMagLength();
        
            % Number of elements (n)
            n = length(x_positionsUI);  % Ensure this matches your model
        
            % Set lower and upper bounds for the magnetic field
            lb = [-0.0025, -0.0025]; % Lower bounds (adjust as needed)
            ub = [0.0025, 0.0025];   % Upper bounds (adjust as needed)
            
            % Set optimization options
            options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
                                   'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 500);
        
            % Define the objective function
            objectiveFunc = @(field) obj.objective_function(field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,obj.world.getOri());
            
            % Perform the optimization
            [optimal_field, optimal_error] = fmincon(objectiveFunc, initial_field, [], [], [], [], lb, ub, [], options);
        
            % If the optimized field has a better fitness, update the best individual
            current_fitness = obj.bestIndividual.getFitness();
            if optimal_error < current_fitness
                obj.bestIndividual.setField(optimal_field);
                disp('Optimized field found with better fitness.');
            else
                disp('No better field found during optimization.');
            end
        end

        function error = objective_function(obj, field, x_positionsUI, y_positionsUI, n, L, magn_momentsUI,orientation)
            % Set the magnetic field
            Bx = field(1);
            By = field(2);
        
            % Solve the elastica using the Solver2D function
            [x_positions, y_positions] = Solver2D(n, L, Bx, By, magn_momentsUI,orientation);
        
            % Calculate the error (Euclidean distance)
            error = sqrt(sum((x_positions - x_positionsUI).^2 + (y_positions - y_positionsUI).^2));
        end

    end



end