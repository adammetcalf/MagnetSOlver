classdef World
    %This class contains the world, including a tentacle, gravity, external
    %magentic fields and the optimisation solver
    
    properties
        Gravity double = [0,0,-9.81];           % Gravity definition
        Tentacle Tentacle;                      % Tentacle object

        % Magnetic Properties
        mu0 double = 4*pi*1e-7;                 % Permeability of free space
        threshold double = 0.001;               % Singularity threshold
        HomegeneousField double = [0,0,0];      % initialise Filed        

        % Magnetic Field Template
        x double;
        y double;
        z double;

        % Magnetic Field
        Bx double;
        By double;
        Bz double;
    end
    
    %% Public Methods
    methods (Access = public)

        %% Constructor
        function obj = World(LinkLength, Angles, Magnetisation, HomogeneousField)
            %WORLD Construct an instance of World

            % Create a tentacle
            obj.Tentacle = Tentacle(LinkLength,Angles,Magnetisation);

            % Create the magnetic field template
            xLow = -0.1;   yLow = -0.1;   zLow = -0.1;
            xHigh = 0.1;   yHigh = 0.1;   zHigh = 0.1;
            xIncr = 0.005;  yIncr = 0.005;  zIncr = 0.005; %steps size

            % Inject homegeneous field in case of re-initialisation
            obj.HomegeneousField = HomogeneousField;
            
            % Create the meshgrid used for magnetic simulation
            [obj.x, obj.y, obj.z] = meshgrid(xLow:xIncr:xHigh, yLow:yIncr:yHigh, zLow:zIncr:zHigh);

            % Init magfield
            obj = obj.InitMagField();

        end
        
        %% Function to update the tentacle joint angles
        function obj = UpDateAngles(obj, Angles)
           
            %   Update the tentacle joint angles
            obj.Tentacle = obj.Tentacle.UpdateAngles(Angles);

            % Re-init magfield
            obj = obj.InitMagField();
        end

        %% Function to Plot
        function obj = plotWorld(obj,PlotOrientationOn,MagField)
            
            %Get the tentacle Joint positions
            HGMs = obj.Tentacle.getHGMs();
            
            for i = 1:size(HGMs,3)
                jointPos(:,i) = HGMs(1:3,4,i);  
            end

            %Get the tentacle Link positions
            LinkPos = obj.Tentacle.getLinks();

            % Get Magnetic moment vectors
            if PlotOrientationOn
                Moments = obj.Tentacle.getMagneticMoments();
            end

            %Check if figure is open, and close if so
            figHandle = findobj('Type', 'figure', 'Number', 1);
                if ~isempty(figHandle)
                    close(figHandle);
                end

            figure(1)
            plot3([jointPos(1,:)],...
                  [jointPos(2,:)],...
                  [jointPos(3,:)],...
                   'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth', 1)
            hold on
            plot3([LinkPos(1,:)],...
                  [LinkPos(2,:)],...
                  [LinkPos(3,:)],...
                   'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b')
            if PlotOrientationOn
                obj.PlotMagenticMoments(LinkPos,Moments)
            end
            if MagField
                quiver3(obj.x, obj.y, obj.z, obj.Bx, obj.By, obj.Bz, 0.5, 'b','LineWidth', 0.25);
            end
            axis equal;
            axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
            grid on
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            hold off
        end

    end

    %% Private Methods
    methods (Access = private)

        %Function to plot the magnetic moments
        function obj = PlotMagenticMoments(obj,LinkPos,Moments)

            for i = 1:size(LinkPos,2)
                Pos = LinkPos(:,i);
                Moment = Moments(:,i);
                quiver3(Pos(1), Pos(2), Pos(3), Moment(1), Moment(2), Moment(3), 'k', 'LineWidth', 1)
            end

        end

        % Function to determine magnetic contribution of each link
        function obj = DetermineMagneticContribution(obj)
                 
            %Get the tentacle Link positions
            LinkPos = obj.Tentacle.getLinks();

            % Get Magnetic moment vectors
            Moments = obj.Tentacle.getMagneticMoments();
            
            % Evaluate Magnetic contribution of each link
            for i = 1:size(LinkPos,2)

                Pos = LinkPos(:,i);
                Moment = Moments(:,i);

                %evaluate magnetic field contribution of dipole
                x1 = obj.x - Pos(1);
                y1 = obj.y - Pos(2);
                z1 = obj.z - Pos(3);
                r1 = sqrt(x1.^2 + y1.^2 + z1.^2);
                rx1 = x1./r1; ry1 = y1./r1; rz1 = z1./r1;
                
                BxDipole = obj.mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*rx1 - Moment(1))./r1.^3;
                ByDipole = obj.mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*ry1 - Moment(2))./r1.^3;
                BzDipole = obj.mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*rz1 - Moment(3))./r1.^3;
            
                % Remove singularities for all dipoles
                BxDipole(r1<obj.threshold) = 0; ByDipole(r1<obj.threshold) = 0; BzDipole(r1<obj.threshold) = 0;

                % Update Overall Magnetic Field
                obj.Bx = obj.Bx + BxDipole;
                obj.By = obj.By + ByDipole;
                obj.Bz = obj.Bz + BzDipole;
            end

        end


        % Function to initialise the Magnetic field
        function obj = InitMagField(obj)

            % Reset Magnetic field to zeros.
            obj.Bx = zeros(size(obj.x));
            obj.By = zeros(size(obj.y));
            obj.Bz = zeros(size(obj.z));

            % Create Magnetic field
            obj.Bx = obj.HomegeneousField(1)*ones(size(obj.x));
            obj.By = obj.HomegeneousField(2)*ones(size(obj.y));
            obj.Bz = obj.HomegeneousField(3)*ones(size(obj.z));

            % Determine Magnetic contribution of each link
            obj = DetermineMagneticContribution(obj);

        end
    end
end

