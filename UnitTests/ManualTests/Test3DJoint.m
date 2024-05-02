clear
clc;

for ii = 0:10:90
    for jj = 0:10:90

    % Define Joints
    joint1 = Joint3D(90,0,0,0,0);
    joint2 = Joint3D(0,ii,jj,0,1);
    joint3 = Joint3D(0,0,0,0,1);
    joint4 = Joint3D(0,0,0,0,1);
    
    % Get frames
    Frame1 = joint1.getFrame();
    Frame2 = joint2.getFrame();
    Frame3 = joint3.getFrame();
    Frame4 = joint3.getFrame();
    
    
    Origin = [1 0 0 0
              0 1 0 0
              0 0 1 0
              0 0 0 1];
    
    %Form HGMs
    HGM = Origin*Frame1;
    
    HGMs(:,:,1) = HGM;
    
    HGM = HGM*Frame2;
    
    HGMs(:,:,2) = HGM;
    
    HGM = HGM*Frame3;
    
    HGMs(:,:,3) = HGM;
    
    HGM = HGM*Frame4;
    
    HGMs(:,:,4) = HGM;
    
    
    % Get joint positions
    for i = 1:size(HGMs,3)
        jointPos(:,i) = HGMs(1:3,4,i);  
    end
    
    figure(1)
    plot3([jointPos(1,:)],...
          [jointPos(2,:)],...
          [jointPos(3,:)],...
           'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
    axis equal;
    axis([-2 2 -2 2 0 4]);
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('90 degrees y')
    view(0,90)
        
    figure(2)
    plot3([jointPos(1,:)],...
          [jointPos(2,:)],...
          [jointPos(3,:)],...
           'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
    axis equal;
    axis([-2 2 -2 2 0 4]);
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('90 degrees y')
    view(-45,30)

    figure(3)
    plot3([jointPos(1,:)],...
          [jointPos(2,:)],...
          [jointPos(3,:)],...
           'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
    axis equal;
    axis([-2 2 -2 2 0 4]);
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('90 degrees y')
    view(-90,0)

    figure(4)
    plot3([jointPos(1,:)],...
          [jointPos(2,:)],...
          [jointPos(3,:)],...
           'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
    axis equal;
    axis([-2 2 -2 2 0 4]);
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('90 degrees y')
    view(0,0)

    pause(0.1)

    end

end