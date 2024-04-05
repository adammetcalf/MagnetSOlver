classdef TestTentacleUnitTests < matlab.unittest.TestCase

    methods(Test)
        % Test constructor
        function constructor(testCase)
            angles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0];
            Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
            tentacle = Tentacle(0.01,angles,12,Magdirection);
            anglesOut = tentacle.getJointAngles();
            testCase.verifyEqual(anglesOut,angles)
        end
        %Test udpdate angles
        function UpdateAngles1(testCase)
            angles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0];
            Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
            tentacle = Tentacle(0.01,angles,12,Magdirection);
            anglesOut = tentacle.getJointAngles();
            angles2 = [90,180;0,90;0,90;0,90;0,90;0,90;0,90];
            tentacle = tentacle.UpdateAngles(angles2);
            anglesOut2 = tentacle.getJointAngles();
            testCase.verifyEqual(angles2,anglesOut2)
        end
        %Test udpdate angles
        function UpdateAngles2(testCase)
            angles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0];
            Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
            tentacle = Tentacle(0.01,angles,12,Magdirection);
            anglesOut = tentacle.getJointAngles();
            angles2 = [90,180;0,90;0,90;0,90;0,90;0,90;0,90];
            tentacle = tentacle.UpdateAngles(angles2);
            anglesOut2 = tentacle.getJointAngles();
            testCase.verifyNotEqual(angles,anglesOut2)
        end
    end
end