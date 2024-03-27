classdef TestIndividualUnitTests < matlab.unittest.TestCase

    methods(Test)
        function constructor(testCase)
            % Define world inputs
            JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
            TentacleMagnetisation = 12;
            LinkLength = 0.01;
            HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
            MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.
            
            % Create World
            world = World(LinkLength, JointAngles, TentacleMagnetisation, HomogeneousField, MultipoleActive);


            %spawn World
            individual = Individual(world);

            angles = individual.getAngles();
            testCase.verifyNotEqual(JointAngles,angles);
        end
        function constructor2(testCase)
            % Define world inputs
            JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
            TentacleMagnetisation = 12;
            LinkLength = 0.01;
            HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
            MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.
            
            % Create World
            world = World(LinkLength, JointAngles, TentacleMagnetisation, HomogeneousField, MultipoleActive);


            %spawn World
            individual = Individual(world);

            angles = individual.getAngles();

            size1 = size(JointAngles);
            size2 = size(angles);


            testCase.verifyEqual(size1,size2);
        end
        function updateAngles(testCase)
            % Define world inputs
            JointAngles = [90,90;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
            TentacleMagnetisation = 12;
            LinkLength = 0.01;
            HomogeneousField = [0,0,25e-3]; % Magnetic field strength in Tesla (25 mT) in +z direction
            MultipoleActive = false; % Include the magnetic effects of the tentacle links on the magentic field.
            
            % Create World
            world = World(LinkLength, JointAngles, TentacleMagnetisation, HomogeneousField, MultipoleActive);


            %spawn World
            individual = Individual(world);

            angles1 = individual.getAngles();

            individual = individual.updateAngles(JointAngles, world);

            angles2 = individual.getAngles();

            testCase.verifyEqual(JointAngles,angles2);
        end
    end
end