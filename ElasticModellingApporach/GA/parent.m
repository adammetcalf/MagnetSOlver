classdef parent
    %PARENT This class contains a parent selected from the population for
    %breeding
    
    properties (Access = private)
        progenator individual;
    end
    
    methods (Access = public)
        %% Constructor
        function obj = parent(popArray)
            %PARENT Construct an instance of this class and return the
            % parent object
            obj = getParent(obj, popArray);
        end
        
        %% Accessors
        function Field = getField(obj)
            Field = obj.progenator.getField();
        end
    end

    methods (Access = private)
        
        %gets a parent from population
        function obj = getParent(obj, popArray)
            %This function gets a parent from the population
            x=rand;
            if x > 0.5 %use Tournament
                obj = tournament(obj, popArray);
                
            else %use Biased Random
                obj = biasedRandom(obj, popArray);
                
            end
        end
        %biased random function
        function obj = biasedRandom(obj, popArray)
            %This function selects an individual from the population using
            %a biased random method
            
            %popArray is already sorted
            bias = rand; 
            if bias>=0.7 % 30% chance
                obj.progenator = popArray(1);
            elseif bias>=0.5 % 20% chance
                obj.progenator = popArray(2);
            elseif bias>=0.35 % 15% chance
                obj.progenator = popArray(3);
            elseif bias>=0.25 % 10% chance
                obj.progenator = popArray(4);
            elseif bias>=0.2 % 5% chance
                obj.progenator = popArray(5);
            elseif bias>=0.15 % 5% chance
                obj.progenator = popArray(6);
            elseif bias>=0.1 % 5% chance
                obj.progenator = popArray(7);
            elseif bias>=0.05 % 5% chance
                obj.progenator = popArray(8);
            elseif bias>=0.025 % 2.5% chance
                obj.progenator = popArray(9);
            else  % 2.5% chance
                obj.progenator = popArray(10);
            end
        end

        %tournament function
        function obj = tournament(obj, popArray)
            %this function randomly selects 2 individuals from the
            %population and selects the fittest one as the progenator
            choice = randi([1 10],1,2);
            while choice(1)==choice(2)
                    choice(2) = randi([1 10],1,1);
            end
            
            %can select fittest from Array position, since array is sorted
            if choice(1)<=choice(2)
                obj.progenator = popArray(choice(1));
            else
                obj.progenator = popArray(choice(2));
            end
            
        end

    end
end
