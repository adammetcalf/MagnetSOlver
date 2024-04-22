close all;
clear;
clc;

% Declare the global variable to retrieve the result
global result;

% Initialize or clear the global variable
result = [];

fig = simpleApp; % This will open the GUI

% Wait for the UI to close
waitfor(fig);

% Check the global variable for the result
%if isempty(result)
%    disp('No result was set.');
%else
%    disp(['Result from GUI: ', result]);
%end

result