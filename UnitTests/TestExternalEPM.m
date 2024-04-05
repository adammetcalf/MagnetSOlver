classdef TestExternalEPM < matlab.unittest.TestCase
    
    
    methods(Test)
       
        function constructor(testCase)
            location = [0;0;1];
            orientation = [1,0,0;0,1,0;0,0,1];
            magnet = ExternalEPM(location,orientation);
            [position,~] = magnet.getMag();
            testCase.verifyEqual(location,position);
        end
        function constructor2(testCase)
            location = [0;0;1];
            orientation = [1,0,0;0,1,0;0,0,1];
            magnet = ExternalEPM(location,orientation);
            [~,moment] = magnet.getMag();
            testCase.verifyEqual(moment,[0;0;970.1]);
        end
    end
    
end