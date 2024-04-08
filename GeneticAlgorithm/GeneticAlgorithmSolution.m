classdef GeneticAlgorithmSolution
    % This class solves the tentacle modelling problem (hopefully)
    
    properties (Access = private)
        world World;
        bestIndividual Individual;
        NextBestIndividual Individual;
    end
    
    methods (Access = public)
        %% Constructor
        function obj = GeneticAlgorithmSolution(world)
            %GENETICALGORITHMSOLUTION Construct an instance of this class

            % init variables
            obj.world = world;
            obj.bestIndividual = Individual(world);
            obj.NextBestIndividual = Individual(world);

        end

        function BestAngles = Train(obj)
                    
            %init epoch count
            epoch = 0;        
            
            % init angle change count
            AngleNoChange = 0;

            while ((epoch < 10) && (AngleNoChange <= 3))
                disp("Epoch: "+ num2str(epoch));

                % Best angles at start of epoch
                AngleStart = obj.bestIndividual.getAngles();

                % perform 1 epoch of training
                obj = obj.performEpoch();

                % Optimise best individual
                obj = obj.Optimise();

                % increment epoch count
                epoch = epoch+1;

                % Best angles at end of epoch
                AngleEnd = obj.bestIndividual.getAngles();

                % Compare best angles at start and end of epoch
                if AngleStart == AngleEnd
                    AngleNoChange = AngleNoChange+1;
                else
                    AngleNoChange =0;
                end

                disp(AngleNoChange);

            end
            
            % Ensure that best angles have propogated through (probably
            % unnecessary)
            BestAngles = obj.bestIndividual.getAngles();

        end
        
    end

    methods (Access = private)
        
        % perform the epoch
        function obj = performEpoch(obj)
            
            % initialise epoch
            pop = Population(obj.world);
            
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

            convCheckCount =0;
            while convCheckCount < 20  
                

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
            % for 20 consequtive evolutions). 

            %update the best individual
            popBest = pop.getBestIndividual();
            popFit = popBest.getFitness();
            oldFit = obj.bestIndividual.getFitness();

            if popFit < oldFit
                obj.bestIndividual = pop.getBestIndividual();
            end

        end


        % Optimisation
        function obj = Optimise(obj)

            % update next best as best. Both will be injected into the
            % population for the next epoch.
            obj.NextBestIndividual = obj.bestIndividual;

            %Get Best individual angles
            angles = obj.bestIndividual.getAngles();

            %Update World angles
            obj.world = obj.world.UpDateAngles(angles);

            %Create PoseSolver
            poseSolver = PoseSolver(obj.world);

            %optimise the angles
            optimizedAngles = poseSolver.optimizeJoints();

            % Set any extreme angles to  zero (unless its
            % Joint0). This prevents 'folding' of the tentacle.
            for i = 1:(size(optimizedAngles,1)-2)

                if (optimizedAngles(i+1,2) >= 179) && (optimizedAngles(i+1,2) <= 181)
                    optimizedAngles(i+1,2) = 0;
                elseif (optimizedAngles(i+1,2) <= -179) && (optimizedAngles(i+1,2) >= -181)
                    optimizedAngles(i+1,2) = 0;   
                else

                end

            end

            % Update Best individual with optimised angles
            obj.bestIndividual = obj.bestIndividual.updateAngles(optimizedAngles,obj.world);

            % Note that the best individual will be automatically injected
            % into the population at the start of the next epoch

        end

    end

end

