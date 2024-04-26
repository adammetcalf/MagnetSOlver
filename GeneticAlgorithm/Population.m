classdef Population
    %POPULATION contains a sorted population of individuals
    
    properties
        popArray Individual; %Array of Individuals
        BestIndividual Individual; %Fittest individual in the population
        AverageFitness double; %Mean fitness of the population
        CurrentMutation;
    end
    
   methods (Access = public)

        %% Constructor
        function obj = Population(World,CurrentMutation)
            obj.CurrentMutation = CurrentMutation;
            for count =1:10
                obj.popArray(count) = Individual(World); %create an individual, and place into the population            
            end
            obj = sortFitness(obj);
        end

        %% Evolve
        function obj = Evolve(obj,World,Evolution)

            %init offspring
            offspringArray = [obj.popArray,obj.popArray];

            % Get 10 parents and perform crosssover
            offspringArray = getParents(obj, offspringArray, World);   

            for i = min(floor((Evolution-1) / 10) + 1, 8)
                % Mutate
                [obj,offspringArray] = Mutate(obj,offspringArray, World, Evolution);
            end

            % Combine popArray and Offsrping Array
            obj.popArray = [obj.popArray,offspringArray];


            % CullArray
            obj = CullArray2(obj,World);


        end

        % function to inject the best individuals from a previous epoch
        function obj = injectBest(obj, individual,NextIndividual)

            %inject the best individuals from previous epoch
            obj.popArray(1) = individual;
            obj.popArray(2) = NextIndividual;

            % Inject individual with reversed angle 1 at end of population
            % (in case of local minima)
            individual = individual.reverseAngleOne();
            obj.popArray(10) = individual;

            obj = sortFitness(obj);

        end


        %% Accessors
        function [Angles, fitness, AverageFitness] = getBest(obj)
            
            %return Avergae fitness of the population
            AverageFitness= obj.AverageFitness;

            %return Angles of the fittest individual
            Angles = obj.BestIndividual.getAngles();

            %return fitness of the fittest individual
            fitness = obj.BestIndividual.getFitness();

        end

        function popArray = getPopArray(obj)
            %return popArray
            popArray = obj.popArray;

        end

        function BestIndividual = getBestIndividual(obj)
            BestIndividual = obj.BestIndividual;

        end
        
    end

   methods (Access = private)

       function offspringArray = getParents(obj, offspringArray,World)

           % get a pair of parents until offspring array is fully replaced
           for count = 1:5

               %Get 2 parents
               mother = Parent(obj.popArray);
               father = Parent(obj.popArray);
               %while mother==father,select father again to ensure no
               %duplication
                while mother.getAngles() == father.getAngles()
                    father = Parent(obj.popArray);
                end

            % extract individual from parent
            mother = mother.getAngles();
            father = father.getAngles();

            % Perfom Crossover
            [obj,child1, child2] = crossOver(obj,mother,father,World);

            % Replace items in offspringarray
            index1 = (2*count)-1;
            index2 = (2*count);
            offSpringArray(index1) =child1;
            offSpringArray(index2) =child2;
            count=count+1;
            end

       end 

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

        function obj = CullArray(obj,World)
            % This function removes the second half of the array

            %must remove duplicates
            pop = obj.popArray;

            index =1; %outer loop iteration count
            for x = 1:length(obj.popArray) %check once for each individual
                ind = pop(index); %get individual
                CompRoute = ind.getAngles(); %get angles for comparison
                
                count = 0; %set duplication count to 0
                innerIndex = 1; %inner loop iteration count
                for y = 1:1:length(obj.popArray) %check population with each individual
                indtoCheck = pop(innerIndex); %get individual to check   
                CheckRoute = indtoCheck.getAngles(); %get angles for comparison
                
                if CompRoute == CheckRoute
                    count = count+1; %if duplicate found, increase count
                    if count > 1 %always expect to find one duplicate, any more replace duplicate
                        pop(innerIndex)= Individual(World); %create an individual 
                    end
                end

                innerIndex = innerIndex + 1;
                end
     
            index = index+1;                    
            end

            obj.popArray = pop;
            newLength= 10;
            
            obj = sortFitness(obj); % order from best to worst
            obj.popArray = obj.popArray(1,1:newLength); % Remove worst half

        end

        % MUCH MORE EFFICIENT CULL.
        function obj = CullArray2(obj, World)
            pop = obj.popArray;
            uniqueAngles = containers.Map('KeyType', 'char', 'ValueType', 'any');
            duplicates = false(1, length(pop)); % Tracks duplicates
        
            for i = 1:length(pop)
                ind = pop(i);
                angleKey = mat2str(ind.getAngles()); % Convert angles to string to use as key
        
                if uniqueAngles.isKey(angleKey)
                    duplicates(i) = true; % Mark as duplicate
                else
                    uniqueAngles(angleKey) = ind; % Add to unique tracker
                end
            end
        
            % Replace duplicates
            for i = find(duplicates)
                pop(i) = Individual(World); % Replace duplicate with new individual
            end
            
            obj.popArray = pop; % Update the object's population array
            newLength= 10;
            obj = sortFitness(obj); % order from best to worst
            obj.popArray = obj.popArray(1,1:newLength); % Remove worst half
        end


        %Crossover Method
        function [obj, child1, child2] = crossOver(obj, mother, father, World)
            % Determine the size of the parent matrices
            n = size(mother, 1);
            
            % Determine how many rows to get from each parent
            choice = randi([floor(n/3), ceil(2*n/3)], 1, 1);
            
            % Create child1 by taking the first 'choice' rows from mother and the remaining rows from father
            child1 = [mother(1:choice, :); father(choice+1:end, :)];
            
            % Ensure we don't exceed the matrix dimensions when choice is equal to n-1
            if choice < n
                % Create child2 by taking the first 'choice' rows from father and the remaining rows from mother
                child2 = [father(1:choice, :); mother(choice+1:end, :)];
            else
                % If choice is n-1, child2 will be exactly like the first 'choice' rows from father followed by the last row from mother
                child2 = [father(1:choice, :); mother(end, :)];
            end

            %Finally, spawn an indivdual for each child using the angles
            temp = Individual(World);
            child1 = temp.updateAngles(child1, World);

            temp = Individual(World);
            child2 = temp.updateAngles(child2, World);

        end

        %Mutation method
        function [obj,offSpringArray] = Mutate(obj,offSpringArray,World, Evolution)
           %this method performs mutation

           selection = randi([1 length(offSpringArray)],1,1); %select random offspring to mutate
           mutant = offSpringArray(selection); %get offspring as individual
           angles = mutant.getAngles(); %get angles
           
           
        % Generate a random integer between 1 and 6
        mutationChoice = randi([1, 6]);
        
        % Choose the mutation method based on the random integer
        switch mutationChoice
            case 1
                angles = obj.mutateSwapRows(angles);
                obj.CurrentMutation.Text = 'Swap Rows';
                drawnow();
            case 2
                angles = obj.mutateRandomResetting(angles);
                obj.CurrentMutation.Text = 'Random Reset';
                drawnow();
            case 3
                angles = obj.mutateCreep(angles);
                obj.CurrentMutation.Text = 'Creep';
                drawnow();
            case 4
                angles = obj.mutateGaussian(angles);
                obj.CurrentMutation.Text = 'Guassian';
                drawnow();
            case 5
                angles = obj.mutateInversion(angles);
                obj.CurrentMutation.Text = 'Inversion';
                drawnow();
            case 6
                angles = obj.mutateUniform(angles, Evolution); % Note: This one uses the Evolution parameter
                obj.CurrentMutation.Text = 'Uniform Mutation';
                drawnow();
        end
          
           mutant = mutant.updateAngles(angles,World); %update angles of individual
           offSpringArray(selection) = mutant;
           
        end

        %% Mutations

        % Scramble mutation
        function angles = mutateSwapRows(obj,angles)

            rLength = size(angles,1); %See how many rows there are
            
            %choose 2 rows in matrix
            choice = randi([1 rLength],1,2); 
            while choice(1)==choice(2)
                 choice(2) = randi([1 rLength],1,1); %ensure the choices are different elements
            end
            
            %perform swap
            y1 = angles(choice(1),:); 
            y2 = angles(choice(2),:);
            angles(choice(1),:) = y2;
            angles(choice(2),:) = y1;

        end

        % Random reset
        function angles = mutateRandomResetting(obj,angles)

            rLength = size(angles,1); %See how many rows there are

            choice = randi([1 rLength],1,1); %select random row

            % Randomly reset the angle of a selected row
            angles(choice, 2) = rand()*360 - 180; % Values between -180 and 180 degrees
        end

        % Creep Mutation
        function angles = mutateCreep(obj, angles)
            % Adds a small random value to the selected angle within a predefined range
            
            rLength = size(angles,1); %See how many rows there are
            choice = randi([1 rLength],1,1); %select random row
            creepRange = 30;

            creepValue = (rand()*2*creepRange) - creepRange; % Creep range defines the max change
            angles(choice, 2) = mod(angles(choice, 2) + creepValue, 360);
            
            if angles(choice, 2) > 180
                angles(choice, 2) = angles(choice, 2) - 360;
            end

        end

        % Gaussian Mutation
        function angles = mutateGaussian(obj, angles)
            % Adds a Gaussian-distributed random value to the selected angle

            stdDev = 60;

            rLength = size(angles,1); %See how many rows there are
            choice = randi([1 rLength],1,1); %select random row

            gaussValue = randn()*stdDev;
            angles(choice, 2) = mod(angles(choice, 2) + gaussValue, 360);
            if angles(choice, 2) > 180
                angles(choice, 2) = angles(choice, 2) - 360;
            end
        end

        % Inversion Mutation
        function angles = mutateInversion(obj, angles)
            % Invert the order of rows within a selected segment

            rLength = size(angles,1); %See how many rows there are
            
            %choose 2 rows in matrix
            choice = randi([1 rLength],1,2); 
            while choice(1)==choice(2)
                 choice(2) = randi([1 rLength],1,1); %ensure the choices are different elements
            end

            if choice(1) < choice(2)
                angles(choice(1):choice(2), :) = flipud(angles(choice(1):choice(2), :));
            end
        end

        % Uniform Mutation
        function angles = mutateUniform(obj, angles, Evolution)

            mutationProb = 0.40; % 40% chance

            % Change the angles in all rows with a certain probability
            for row = 1:size(angles, 1)
                if rand() < mutationProb
                    angles(row, 2) = rand()*360 - 180; % Values between -180 and 180 degrees
                end
            end
        end

    end
end

