classdef population
    % population contains a sorted population of individuals

    %% private properties
    properties (Access = private)
        popArray individual; % array of individuals
        BestIndividual individual; % Fittest individual in the population
        AverageFitness double; % Mean fitness of the population
    end

    %% public methods
    methods (Access = public)

        %% Constructor
        function obj = population(world)
            
            % create an array of individuals
            for count = 1:20
                obj.popArray = individual(world);
            end
            obj = sortFitness(obj);

        end

        %% Evolve
        function obj = Evolve(obj,world,evolution)

            %init offspring
            offspringarray = [obj.popArray,obj.popArray];

            % Get 10 parents and perform corssover
            offspringarray = getParents(obj, offspringarray, world);

            for i = min(floor((evolution-1) / 10) + 1, 8)
                % Mutate
                [obj,offspringArray] = Mutate(obj,offspringArray, world);
            end

            % Combine popArray and Offsrping Array
            obj.popArray = [obj.popArray,offspringArray];


            % CullArray
            obj = CullArray2(obj, world);

        end



        %% Accessors
        function popArray = getPopArray(obj)
            
            popArray = obj.popArray;
        end

        function BestIndividual = getBestIndividual(obj)

            BestIndividual = obj.BestIndividual;
        end

        function [Field, fitness, AverageFitness] = getBest(obj)
            
            %return Avergae fitness of the population
            AverageFitness= obj.AverageFitness;

            %return the field of the fittest individual
            Field = obj.BestIndividual.getField();

            %return fitness of the fittest individual
            fitness = obj.BestIndividual.getFitness();

        end

    end

    %% private methods
    methods (Access = private)

        function obj = sortFitness(obj)
            %This method sorts the population array from best to worst
            %fitness
            
            sizeSort = length(obj.popArray); %find how many individuals in pop

            %for as many individuals are in the population, create a
            %temporary array of just the fitness scores
            for i = 1:sizeSort
                tempIndividual = obj.popArray(i);
                fitnessArray(i) = tempIndividual.getFitness();
            end

            %sort fitnessArray (B is the array sorted lowest to highest, I
            %is the array index from fitness array)
            [B,I] = sort(fitnessArray);
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


        function offspringArray = getParents(obj, offspringArray,world)

           % get a pair of parents until offspring array is fully replaced
           for count = 1:10

               %Get 2 parents
               mother = parent(obj.popArray);
               father = parent(obj.popArray);

                   %while mother==father,select father again to ensure no
                   %duplication
                    while mother.getField() == father.getField()
                        father = parent(obj.popArray);
                    end

                % extract individual from parent
                mother = mother.getField();
                father = father.getField();
    
                % Perfom Crossover
                [obj,child1, child2] = crossOver(obj,mother,father,world);
    
                % Replace items in offspringarray
                index1 = (2*count)-1;
                index2 = (2*count);
                offSpringArray(index1) =child1;
                offSpringArray(index2) =child2;
                count=count+1;
            end

        end

        %Crossover Method
        function [obj, child1, child2] = crossOver(obj, mother, father, world)


            % Create children by combining fields
            child1 = [mother(1); father(2)];
            child2 = [father(1); mother(2)];
            
            %Finally, spawn an indivdual for each child using the angles
            temp = individual(World);
            child1 = temp.updateField(child1, world);

            temp = Individual(World);
            child2 = temp.updateField(child2, world);

        end

        %Mutation method
        function [obj,offSpringArray] = Mutate(obj,offSpringArray, world)
           %this method performs mutation

           selection = randi([1 length(offSpringArray)],1,1); %select random offspring to mutate
           mutant = offSpringArray(selection); %get offspring as individual
           Field = mutant.getField(); %get Field
           
           
        % Generate a random integer between 1 and 6
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
          
           mutant = mutant.updateField(Field,world); %update angles of individual
           offSpringArray(selection) = mutant;
           
        end

        % MUCH MORE EFFICIENT CULL.
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
            newLength= 20;
            obj = sortFitness(obj); % order from best to worst
            obj.popArray = obj.popArray(1,1:newLength); % Remove worst half
        end

        %% Mutations

        % Scramble mutation
        function Field = mutateSwapRows(obj,Field)

            %Swap Field values
            
            temp = Field;
            Field(1) = temp(1);
            Field(2) = temp(2);

        end

        % Random reset
        function Field = mutateRandomResetting(obj,Field)

            % Randomly reset one of the Fields

            % Randomly select one of the two elements (either index 1 or 2)
            index = randi(2);

            % Define limits
            lower_limit = -0.0025;  
            upper_limit = 0.0025; 
            
            % Randomize a and b between the limits
            a = lower_limit + (upper_limit - lower_limit) * rand;

            Field(index) = a;

        end

        % Creep Mutation
        function Field = mutateCreep(obj, Field)
            % Adds a small random value to the selected Field within a predefined range
            
            % Randomly select one of the two elements (either index 1 or 2)
            index = randi(2);

            % define the range of the creep
            creepRange = 0.0001;

            % find the creep value
            creepValue = (rand()*2*creepRange) - creepRange; % Creep range defines the max change
            
            Field(index) = Field(index)+creepvalue;
            
            % ensure that both elements are still in the valid range
            Field = max(min(Field, 0.0025), -0.0025);
        end

        % Gaussian Mutation
        function Field = mutateGaussian(obj, Field)
            % Adds a Gaussian-distributed random value to the selected
            % Field

            stdDev = 0.0005;

            % Randomly select one of the two elements (either index 1 or 2)
            index = randi(2);

            Field(index) = Field(index)+stdDev;
            
            % ensure that both elements are still in the valid range
            Field = max(min(Field, 0.0025), -0.0025);
        end

        % Inversion Mutation
        function Field = mutateInversion(obj, Field)
            % Invert the Field of selected direction

            % Randomly select one of the two elements (either index 1 or 2)
            index = randi(2);

            Field(index) = -(Field(index));
        end


    end


end
