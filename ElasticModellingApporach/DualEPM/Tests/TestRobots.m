close all;
clear;
clc;

% Define robots
robot1 = RobotWithMagnet('urdf/kuka_iiwa_1.urdf');
robot2 = RobotWithMagnet('urdf/kuka_iiwa_2.urdf');


% Define meshgrid over which to compute and visualize the magnetic field
[x_grid, y_grid, z_grid] = meshgrid(linspace(-0.25, 0.25, 25), ...
                                    linspace(-0.25, 0.25, 25), ...
                                    linspace(-0.25, 0.25, 25));


% Initialize magnetic field components
Bx = zeros(size(x_grid));
By = zeros(size(y_grid));
Bz = zeros(size(z_grid));


[position1, moment1] = robot1.getMagneticField();
[position2, moment2] = robot2.getMagneticField();

% Compute the magnetic field at each point in the grid
for i = 1:numel(x_grid)
    point = [x_grid(i), y_grid(i), z_grid(i)];
    
    % Magnetic field contribution from robot1
    B1 = calculateMagneticField(point, position1, moment1);
    
    % Magnetic field contribution from robot2
    B2 = calculateMagneticField(point, position2, moment2);
    
    % Total magnetic field
    B_total = B1 + B2;
    
    % Store the field components
    Bx(i) = B_total(1);
    By(i) = B_total(2);
    Bz(i) = B_total(3);
end

% Visualization of robots and tentacle
figure(1);
hold on;
robot1.showRobot();
robot2.showRobot();
quiver3(x_grid, y_grid, z_grid, Bx, By, Bz);


% Update joint angles
Jointangles =  [-40,0,-40,90,30,-30,-100];

robot1 = robot1.updateConfiguration(Jointangles);
robot2 = robot2.updateConfiguration(Jointangles);

test=robot1.getCurrentConfig();

[position1, moment1] = robot1.getMagneticField();
[position2, moment2] = robot2.getMagneticField();

% Compute the magnetic field at each point in the grid
for i = 1:numel(x_grid)
    point = [x_grid(i), y_grid(i), z_grid(i)];
    
    % Magnetic field contribution from robot1
    B1 = calculateMagneticField(point, position1, moment1);
    
    % Magnetic field contribution from robot2
    B2 = calculateMagneticField(point, position2, moment2);
    
    % Total magnetic field
    B_total = B1 + B2;
    
    % Store the field components
    Bx(i) = B_total(1);
    By(i) = B_total(2);
    Bz(i) = B_total(3);
end

% Visualization of robots and tentacle
figure(2);
hold on;
robot1.showRobot();
robot2.showRobot();
quiver3(x_grid, y_grid, z_grid, Bx, By, Bz);



function B_total = calculateMagneticField(point, position, moment)
    % calculateMagneticField - Computes the magnetic field at a given point 
    % due to a magnetic dipole.
    %
    % Syntax: B_total = calculateMagneticField(point, position, moment)
    %
    % Inputs:
    %    point - The point in space [x, y, z] where the field is calculated
    %    position - The position of the magnetic dipole [x, y, z]
    %    moment - The magnetic moment of the dipole [mx, my, mz]
    %
    % Outputs:
    %    B_total - The magnetic field vector [Bx, By, Bz] at the given point

    % Compute the displacement vector R from the dipole to the point
    R = point(:) - position(:); % Ensure R is a column vector
    
    % Calculate the norm (magnitude) of R
    R_norm = norm(R);
    
    % Calculate the dot product of the magnetic moment and R
    m_dot_R = dot(moment, R);
    
    % Calculate the magnetic field using the dipole formula
    if R_norm ~= 0
        B_total = (3 * R * m_dot_R / R_norm^5) - (moment / R_norm^3);
    else
        B_total = [0, 0, 0]; % Handle the case where point coincides with position
    end
end