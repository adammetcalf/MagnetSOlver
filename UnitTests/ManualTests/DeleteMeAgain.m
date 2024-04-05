close all;
clear;
clc;

% This test script ensures that having a large and powerful external magent
% outside of the bounds of the simulated magnetic workspace can still be
% included in the field.

% Define EPM
location = [0;0;0.15];
orientation = [1,0,0;0,1,0;0,0,1];
magnet = ExternalEPM(location,orientation);
[Pos,Moment] = magnet.getMag();
mu0 = 4*pi*1e-7; 


% Create the magnetic field template
xLow = -0.1;   yLow = -0.1;   zLow = -0.1;
xHigh = 0.1;   yHigh = 0.1;   zHigh = 0.1;
xIncr = 0.005;  yIncr = 0.005;  zIncr = 0.005; %steps size

% Magnetic field strength in Tesla (25 mT) in +z direction
HomogeneousField = [0,0,0]; 

% Create the meshgrid used for magnetic simulation
[x, y, z] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh, zLow:zIncr:zHigh);


% Reset Magnetic field to zeros.
Bx = zeros(size(x));
By = zeros(size(y));
Bz = zeros(size(z));

% Create Magnetic field
Bx = HomogeneousField(1)*ones(size(x));
By = HomogeneousField(2)*ones(size(y));
Bz = HomogeneousField(3)*ones(size(z));

%evaluate magnetic field contribution of dipole
[BxDipole,ByDipole,BzDipole] = calcFieldContribution(Pos,Moment,x,y,z,mu0);

% Update Overall Magnetic Field
Bx = Bx + BxDipole;
By = By + ByDipole;
Bz = Bz + BzDipole;

quiver3(x, y, z, Bx, By, Bz, 0.5, 'b','LineWidth', 0.25);
axis equal;
axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off

% Calculate the FieldContirbution of an individual Dipole
function [BxDipole,ByDipole,BzDipole] = calcFieldContribution(Pos,Moment,x,y,z,mu0)
threshold = 0.0001;
    x1 = x - Pos(1);
    y1 = y - Pos(2);
    z1 = z - Pos(3);
    r1 = sqrt(x1.^2 + y1.^2 + z1.^2);
    rx1 = x1./r1; ry1 = y1./r1; rz1 = z1./r1;
    
    BxDipole = mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*rx1 - Moment(1))./r1.^3;
    ByDipole = mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*ry1 - Moment(2))./r1.^3;
    BzDipole = mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*rz1 - Moment(3))./r1.^3;

    % Remove singularities for dipole
    BxDipole(r1<threshold) = 0; ByDipole(r1<threshold) = 0; BzDipole(r1<threshold) = 0;
end