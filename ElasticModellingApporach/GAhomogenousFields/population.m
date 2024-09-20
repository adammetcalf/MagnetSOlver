classdef population
    % population contains a sorted population of individuals

    %% private properties
    properties (Access = public)
        popArray individual; % array of individuals
        BestIndividual individual; % Fittest individual in the population
        AverageFitness double; % Mean fitness of the population
    end

    %% public methods
    methods (Access = public)

        %% Constructor
        function obj = population(world)
            
            % create an array of individuals
            for count = 1:10
                obj.popArray(count) = individual(world);
            end
            obj = sortFitness(obj);

        end

        %% Evolve
        function obj = Evolve(obj, world, evolution)

            %init offspring
            offspringarray = [obj.popArray, obj.popArray];

            % Get 10 parents and perform crossover
            offspringarray = getParents(obj, offspringarray, world);

            for i = min(floor((evolution-1) / 10) + 1, 8)
                % Mutate
                [obj, offspringarray] = Mutate(obj, offspringarray, world);
            end

            % Combine popArray and Offspring Array
            obj.popArray = [obj.popArray, offspringarray];

            % CullArray
            obj = CullArray2(obj, world);

        end

        %% Inject best
        % function to inject the best individuals from a previous epoch
        function obj = injectBest(obj, individual, NextIndividual)

            %inject the best individuals from previous epoch
            obj.popArray(1) = individual;
            obj.popArray(2) = NextIndividual;

            obj = sortFitness(obj);

        end

        %% Accessors
        function popArray = getPopArray(obj)
            popArray = obj.popArray;
        end

        function BestIndividual = getBestIndividual(obj)
            BestIndividual = obj.BestIndividual;
        end

        function [Field, fitness, AverageFitness] = getBest(obj)
            %return Average fitness of the population
            AverageFitness = obj.AverageFitness;

            %return the field of the fittest individual
            Field = obj.BestIndividual.getField();

            %return fitness of the fittest individual
            fitness = obj.BestIndividual.getFitness();
        end

    end

    %% private methods
    methods (Access = private)

        function obj = sortFitness(obj)
            %This method sorts the population array from best to worst fitness
            sizeSort = length(obj.popArray); %find how many individuals in pop

            %for as many individuals are in the population, create a
            %temporary array of just the fitness scores
            for i = 1:sizeSort
                tempIndividual = obj.popArray(i);
                fitnessArray(i) = tempIndividual.getFitness();
            end

            %sort fitnessArray (B is the array sorted lowest to highest, I
            %is the array index from fitness array)
            [B, I] = sort(fitnessArray);
            obj.AverageFitness = mean(B);
            
            %use newly found I to build temp array of ordered individuals
            for i = 1:sizeSort
                fromPopIndex = I(i);
                tempArray(i) = obj.popArray(fromPopIndex);
            end
            
            % Update population array 
            obj.popArray = tempArray;

            % update best individual
            obj.BestIndividual = obj.popArray(1);
        end

        function offspringArray = getParents(obj, offspringArray, world)
            % get a pair of parents until offspring array is fully replaced
            for count = 1:5
                % Get 2 parents
                mother = parent(obj.popArray);
                father = parent(obj.popArray);

                % while mother == father, select father again to ensure no duplication
                while mother.getField() == father.getField()
                    father = parent(obj.popArray);
                end

                % extract individual from parent
                mother = mother.getField();
                father = father.getField();

                % Perform Crossover
                [obj, child1, child2] = crossOver(obj, mother, father, world);

                % Replace items in offspringarray
                index1 = (2 * count) - 1;
                index2 = (2 * count);
                offspringArray(index1) = child1;
                offspringArray(index2) = child2;
                count = count + 1;
            end
        end

        %Crossover Method
        function [obj, child1, child2] = crossOver(obj, mother, father, world)
            
            % Create children by combining fields
            alpha = rand;
            child1 = alpha * mother + (1 - alpha) * father;
            child2 = (1 - alpha) * mother + alpha * father;
            
            %Finally, spawn an individual for each child using the angles
            temp = individual(world);
            child1 = temp.updateField(child1, world);

            temp = individual(world);
            child2 = temp.updateField(child2, world);
        end

        %Mutation method
        function [obj, offspringArray] = Mutate(obj, offspringArray, world)
            %this method performs mutation
            selection = randi([1 length(offspringArray)], 1, 1); %select random offspring to mutate
            mutant = offspringArray(selection); %get offspring as individual
            Field = mutant.getField(); %get Field
            
            % Generate a random integer between 1 and 5
            mutationChoice = randi([1, 5]);
        
            % Choose the mutation method based on the random integer
            switch mutationChoice
                case 1
                    Field = obj.mutateSwapRows(Field);
                case 2
                    Field = obj.mutateRandomResetting(Field);
                case 3
                    Field = obj.mutateCreep(Field);
                case 4
                    Field = obj.mutateGaussian(Field);
                case 5
                    Field = obj.mutateInversion(Field);
            end
          
            mutant = mutant.updateField(Field, world); %update angles of individual
            offspringArray(selection) = mutant;
        end

        % Much more efficient cull.
        function obj = CullArray2(obj, world)
            pop = obj.popArray;
            uniqueFields = containers.Map('KeyType', 'char', 'ValueType', 'any');
            duplicates = false(1, length(pop)); % Tracks duplicates
        
            for i = 1:length(pop)
                ind = pop(i);
                fieldKey = mat2str(ind.getField()); % Convert Fields to string to use as key
        
                if uniqueFields.isKey(fieldKey)
                    duplicates(i) = true; % Mark as duplicate
                else
                    uniqueFields(fieldKey) = ind; % Add to unique tracker
                end
            end
        
            % Replace duplicates
            for i = find(duplicates)
                pop(i) = individual(world); % Replace duplicate with new individual
            end
            
            obj.popArray = pop; % Update the object's population array
            newLength = 20;
            obj = sortFitness(obj); % order from best to worst
            obj.popArray = obj.popArray(1:newLength); % Remove worst half
        end

        %% Mutations

        % Scramble mutation
        function Field = mutateSwapRows(obj, Field)
            % Swap Field values
            temp = Field;
            Field(1) = temp(1);
            Field(2) = temp(2);
        end

        % Random reset
        function Field = mutateRandomResetting(obj, Field)
            % Randomly reset one of the Fields
            index = randi(2); % Randomly select one of the two elements
            lower_limit = -0.0025;  
            upper_limit = 0.0025; 
            a = lower_limit + (upper_limit - lower_limit) * rand;
            Field(index) = a;
        end

        % Creep Mutation
        function Field = mutateCreep(obj, Field)
            % Adds a small random value to the selected Field within a predefined range
            index = randi(2); % Randomly select one of the two elements
            creepRange = 0.0001;
            creepValue = (rand() * 2 * creepRange) - creepRange; % Creep range defines the max change
            Field(index) = Field(index) + creepValue;
            Field = max(min(Field, 0.0025), -0.0025); % ensure that both elements are still in the valid range
        end

        % Gaussian Mutation
        function Field = mutateGaussian(obj, Field)
            % Adds a Gaussian-distributed random value to the selected Field
            stdDev = 0.0005;
            index = randi(2); % Randomly select one of the two elements
            Field(index) = Field(index) + stdDev;
            Field = max(min(Field, 0.0025), -0.0025); % ensure that both elements are still in the valid range
        end

        % Inversion Mutation
        function Field = mutateInversion(obj, Field)
            % Invert the Field of selected direction
            index = randi(2); % Randomly select one of the two elements
            Field(index) = -(Field(index));
        end

    end

end

