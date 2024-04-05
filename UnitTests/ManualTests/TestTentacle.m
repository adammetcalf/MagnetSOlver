close all;
clear;
clc;


angles = [90,180;0,0;0,0;0,0;0,0;0,0;0,0];
Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];

Tentacle = Tentacle(0.01,angles,12,Magdirection);

HGMs = Tentacle.getHGMs();

J1 = Tentacle.getJacobean();

Test = Tentacle.getLinkHGMs();

for i = 1:size(HGMs,3)
    Position(:,i) = HGMs(1:3,4,i);  
end

figure(1)
plot3([Position(1,:)],...
      [Position(2,:)],...
      [Position(3,:)],...
       'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
hold on
axis equal;
axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')


angles = [90,45;0,0;0,0;0,0;0,0;0,0;0,0];

Tentacle = Tentacle.UpdateAngles(angles);

HGMs = Tentacle.getHGMs();

J2 = Tentacle.getJacobean();

for i = 1:size(HGMs,3)
    Position(:,i) = HGMs(1:3,4,i);  
end

figure(2)
plot3([Position(1,:)],...
      [Position(2,:)],...
      [Position(3,:)],...
       'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
hold on
axis equal;
axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')