close all;
clear;
clc;

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.TestReportPlugin

suiteFolder = TestSuite.fromFolder("UnitTests");
result = run(suiteFolder);

numPassed = sum([result.Passed]);

% Number of failed tests
numFailed = sum([result.Failed]);

% Display the results
fprintf('Number of Passed Tests: %d\n', numPassed);
fprintf('Number of Failed Tests: %d\n', numFailed);