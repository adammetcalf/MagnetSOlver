classdef TestPopulationUnitTests < matlab.unittest.TestCase

    methods(Test)
        function constructor(testCase)
            % Define world inputs
            JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
            Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
            TentacleMagnetisation = 12;
            LinkLength = 0.01;
            HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
            MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.
            
            % Create World
            world = World(LinkLength, JointAngles, TentacleMagnetisation, Magdirection, HomogeneousField,[], MultipoleActive);

            %Create population
            pop = Population(world);

            %get population array
            popArray = pop.getPopArray();

            l = length(popArray);

            testCase.verifyEqual(l,10);
        end
        % Check poparray is sorted in fitness order (smallest fitness ==
        % smaller index)
        function constructor2(testCase)
            % Define world inputs
            JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
            Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
            TentacleMagnetisation = 12;
            LinkLength = 0.01;
            HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
            MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.
            
            % Create World
            world = World(LinkLength, JointAngles, TentacleMagnetisation, Magdirection, HomogeneousField,[], MultipoleActive);

            %Create population
            pop = Population(world);

            %get population array
            popArray = pop.getPopArray();

            l = length(popArray);

            TestAns = [false,false,false,false,false,false,false,false,false];

            for i=1:(l-1)
                ind1 = popArray(i);
                ind2 = popArray(i+1);

                fit1 = ind1.getFitness();
                fit2 = ind2.getFitness();

                if (fit1 <= fit2)
                    TestAns(i) = true;
                else
                    TestAns(i) = false;
                end

            end

            DesiredAns = [true,true,true,true,true,true,true,true,true];

            testCase.verifyEqual(l,10);
        end
        function Evolve(testCase)
            % Define world inputs
            JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
            Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
            TentacleMagnetisation = 12;
            LinkLength = 0.01;
            HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
            MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.
            
            % Create World
            world = World(LinkLength, JointAngles, TentacleMagnetisation, Magdirection, HomogeneousField,[], MultipoleActive);

            %Create population
            pop = Population(world);

            %get population array
            popArray = pop.getPopArray();

            for i=1:length(popArray)
                ind = popArray(i);
                fitness1(i) = ind.getFitness();

            end

            pop = pop.Evolve(world,25);

            %get population array
            popArray = pop.getPopArray();

            for i=1:length(popArray)
                ind = popArray(i);
                fitness2(i) = ind.getFitness();

            end

            testCase.verifyNotEqual(fitness1,fitness2);
        end

    end
end