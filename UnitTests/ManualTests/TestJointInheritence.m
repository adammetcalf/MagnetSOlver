close all;
clear;
clc;


%Declare empty array of Inetrface

Joints = cell(1);

%try to replace Joint(1) with 2d joint child
joint = Joint(90,180,0,0);
Joints{1} = joint;

%Try to replace with 3d Joint3D child
joint = Joint3D(90,180,0,0,0);
Joints{2} = joint;