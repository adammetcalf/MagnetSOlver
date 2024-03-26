classdef TestJointUnitTests < matlab.unittest.TestCase

    methods(Test)
        function constructor(testCase)
            joint = Joint(0,90,0,1);
            Frame = joint.getFrame();
            expSolution = [1 0 0 0;
                           0 0 -1 0
                           0 1 0 1
                           0 0 0 1];
            testCase.verifyEqual(Frame,expSolution);
        end
        function updater(testCase)
            joint = Joint(0,90,0,1);
            Frame1 = joint.getFrame();
            joint = joint.UpdateAngles([90,90]);
            Frame2 = joint.getFrame();
            testCase.verifyNotEqual(Frame1,Frame2);
        end
        function updater2(testCase)
            joint = Joint(0,90,0,1);
            joint = joint.UpdateAngles([90,90]);
            Frame = joint.getFrame();
            expSolution = [0 0 1 0;
                           1 0 0 0
                           0 1 0 1
                           0 0 0 1];
            testCase.verifyEqual(Frame,expSolution);
        end
    end
end