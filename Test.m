close all;
clear;
clc;


angles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0];

Tentacle = Tentacle(2,angles);

HGMs = Tentacle.HGMs;

for i = 1:size(HGMs,3)
    Position(:,i) = HGMs(1:3,4,i);  
end

