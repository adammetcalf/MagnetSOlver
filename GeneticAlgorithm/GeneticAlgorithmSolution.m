classdef GeneticAlgorithmSolution
    % This class solves the tentacle modelling problem (hopefully)
    
    properties (Access = private)
        world World;
        bestIndividual Individual;
    end
    
    methods (Access = public)
        %% Constructor
        function obj = GeneticAlgorithmSolution(world)
            %GENETICALGORITHMSOLUTION Construct an instance of this class

            % init variables
            obj.world = world;
            obj.bestIndividual = Individual(world);

        end

        function [obj, BestAngles] = Train(obj)

            epoch = 0;
           
            while epoch < 30
                disp("Epoch: "+ num2str(epoch));


                % perform 1 epoch of training
                obj = obj.performEpoch();

                % increment epoch count
                epoch = epoch+1;
            end
            
            BestAngles = obj.bestIndividual.getAngles();


        end
        
    end

    methods (Access = private)
        
        % perform the epoch
        function obj = performEpoch(obj)
            
            % initialise epoch
            pop = Population(obj.world);
            
            %inject best individual from previous epoch
            pop = pop.injectBest(obj.bestIndividual);

            % Perform the evolution cycles
            obj = obj.performEvolution(pop);
 
        end

        % perform the evolution until there is no further convolution
        function obj = performEvolution(obj, pop)

            % init convolution check count, which is used to track
            % improvements across evloution cycles.

            Evolution = 0;

            convCheckCount =0;
            while convCheckCount < 50  
                

                % Get the best individual
                [~, fitness, ~] = pop.getBest();

                % Evolve the population
                pop = pop.Evolve(obj.world, Evolution);

                % Get the best individual after evolution
                [~, fitness2, ~] = pop.getBest();

                % if there has been an improvement after evolution, reset
                % convolution count to zero, else increment
                if fitness2 == fitness
                    convCheckCount = convCheckCount+1;
                else
                    convCheckCount = 0;
                end
                Evolution = Evolution+1;

                disp("Evolution: "+ num2str(Evolution)+ ", Fitness: "+ num2str(fitness2));

            end 

            % End of the evolution cycles for this epoch (ie no improvement
            % for 50 consequtive evolutions). 

            %update the best individual
            popBest = pop.getBestIndividual();
            popFit = popBest.getFitness();
            oldFit = obj.bestIndividual.getFitness();

            if popFit < oldFit
                obj.bestIndividual = pop.getBestIndividual();
            end

        end

    end

end

