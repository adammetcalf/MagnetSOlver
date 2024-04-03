classdef World
    %This class contains the world, including a tentacle, gravity, external
    %magentic fields and the optimisation solver
    
    properties (Access = public)
        Gravity double = -9.81;                 % Gravity definition
        Tentacle Tentacle;                      % Tentacle object
        MagForceTorque double;                  % Matrix to contain the Magnetic forces and torques on each link
        Fg double;                              % Gravitational force on each link
        StiffnessEffects double;                % Matrix to contain the effects of stiffness on each joint
        PlotLength double;                      % Axis limits, related to tentacle length

        % Magnetic Properties
        mu0 double = 4*pi*1e-7;                 % Permeability of free space
        threshold double = 0.001;               % Singularity threshold
        HomegeneousField double = [0,0,0];      % initialise Field  
        IncludeMultiPole logical = false;       % Include effects of link magnetism in overall field

        % Magnetic Field Template
        x double;
        y double;
        z double;

        % Magnetic field discrete step sizes
        xIncr double = 0;
        yIncr double = 0;
        zIncr double = 0;

        % Magnetic Field
        Bx double;
        By double;
        Bz double;

        
    end
    
    %% Public Methods
    methods (Access = public)

        %% Constructor
        function obj = World(LinkLength, Angles, Magnetisation, HomogeneousField,MultipoleActive)
            %WORLD Construct an instance of World

            % Create a tentacle
            obj.Tentacle = Tentacle(LinkLength,Angles,Magnetisation);

            % Get overall tentacle length
            obj.PlotLength = (size(Angles,1)+2)*LinkLength;

            % Create the magnetic field template
            xLow = -0.1;   yLow = -0.1;   zLow = -0.1;
            xHigh = 0.1;   yHigh = 0.1;   zHigh = 0.1;
            obj.xIncr = 0.005;  obj.yIncr = 0.005;  obj.zIncr = 0.005; %steps size

            % Inject homegeneous field in case of re-initialisation
            obj.HomegeneousField = HomogeneousField;

            obj.IncludeMultiPole = MultipoleActive;
            
            % Create the meshgrid used for magnetic simulation
            [obj.x, obj.y, obj.z] = meshgrid(xLow:obj.xIncr:xHigh, yLow:obj.yIncr:yHigh, zLow:obj.zIncr:zHigh);

            % Init magfield
            obj = obj.InitMagField();

            % Evaluate Magnetic Forces and Torques
            obj = EvaluateMagneticForces(obj);

            % Evaluate Gravtiational Forces
            obj = EvaluateGravity(obj);

        end
        
        %% Function to update the tentacle joint angles
        function obj = UpDateAngles(obj, Angles)
           
            %   Update the tentacle joint angles
            obj.Tentacle = obj.Tentacle.UpdateAngles(Angles);

            % Re-init magfield
            obj = obj.InitMagField();

            % Evaluate Magnetic Forces and Torques
            obj = EvaluateMagneticForces(obj);

            % Evaluate Gravtiational Forces
            obj = EvaluateGravity(obj);
        end

        %% Function to Plot
        function obj = plotWorld(obj,PlotOrientationOn,MagField,figureNumber)
            
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
            figHandle = findobj('Type', 'figure', 'Number', figureNumber);
                if ~isempty(figHandle)
                    close(figHandle);
                end

            figure(figureNumber)
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
            axis([-obj.PlotLength obj.PlotLength -obj.PlotLength obj.PlotLength -obj.PlotLength obj.PlotLength]);
            grid on
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            hold off
        end


        %% Accessors
        %get the joint angles
        function Angles = getJointAngles(obj)
            Angles = obj.Tentacle.getJointAngles();

        end

        % Get the torques and forces
        function ForceTorques = getForcesTorques(obj)

            %MagForceTorques = [Fx;Fy;Fz;Tx;Ty;Tz];

            ForceTorques = obj.MagForceTorque; %  Magnetic Force Torque Matrix
            ForceTorques(3,:) = ForceTorques(2,:) + obj.Fg; % magnetic force torque matrix + effects of gravity acting in Z direction

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
                [BxDipole,ByDipole,BzDipole] = calcFieldContribution(obj,Pos,Moment);

                % Update Overall Magnetic Field
                obj.Bx = obj.Bx + BxDipole;
                obj.By = obj.By + ByDipole;
                obj.Bz = obj.Bz + BzDipole;
            end

        end

        % Calculate the FieldContirbution of an individual Dipole
        function [BxDipole,ByDipole,BzDipole] = calcFieldContribution(obj,Pos,Moment)
                x1 = obj.x - Pos(1);
                y1 = obj.y - Pos(2);
                z1 = obj.z - Pos(3);
                r1 = sqrt(x1.^2 + y1.^2 + z1.^2);
                rx1 = x1./r1; ry1 = y1./r1; rz1 = z1./r1;
                
                BxDipole = obj.mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*rx1 - Moment(1))./r1.^3;
                ByDipole = obj.mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*ry1 - Moment(2))./r1.^3;
                BzDipole = obj.mu0/(4*pi) * (3*(Moment(1)*rx1 + Moment(2)*ry1 + Moment(3)*rz1).*rz1 - Moment(3))./r1.^3;
            
                % Remove singularities for dipole
                BxDipole(r1<obj.threshold) = 0; ByDipole(r1<obj.threshold) = 0; BzDipole(r1<obj.threshold) = 0;
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

            if obj.IncludeMultiPole
                % Determine Magnetic contribution of each link
                obj = DetermineMagneticContribution(obj);
            end

        end

        % Function to evaluate magentic forces at each joint
        function obj = EvaluateMagneticForces(obj)

            %Get the tentacle Link positions
            LinkPos = obj.Tentacle.getLinks();

            % Get Magnetic moment vectors
            Moments = obj.Tentacle.getMagneticMoments();


            % Evaluate Torque/Force at each point
            for i = 1:size(LinkPos,2) 
                
                %Create a local copy of the magnetic field in order to
                %remove the contribution of this Link
                BxLocal = obj.Bx;
                ByLocal = obj.By;
                BzLocal = obj.Bz;

                %Get the link moment and position
                Pos = LinkPos(:,i);
                Moment = Moments(:,i);

                if obj.IncludeMultiPole
                    % Remove the effects of this particular point of interest
                    % from the magnetic field to remove self-interaction
                    % effects
                    [BxDipole,ByDipole,BzDipole] = calcFieldContribution(obj,Pos,Moment);
                    BxLocal = BxLocal-BxDipole;
                    ByLocal = ByLocal-ByDipole;
                    BzLocal = BzLocal-BzDipole;
                    % The Magnetic field contained by BxLocal etc is not
                    % appropriate to apply to this Link
                end

                %Find closest grid point to POI - this is where we will be
                % actually applying the field
                [idx, idy, idz] = findClosestGridPoint(obj, Pos);

                % Compute spacial derivatives
                [dBx_dx, dBx_dy, dBx_dz] = gradient(BxLocal, obj.x(1,:,1), obj.y(:,1,1), obj.z(1,1,:));
                [~, dBy_dy, dBy_dz] = gradient(ByLocal, obj.x(1,:,1), obj.y(:,1,1), obj.z(1,1,:));
                

                % Construct U vector
                U = [dBx_dx(idx, idy, idz); ...
                     dBx_dy(idx, idy, idz); ...
                     dBx_dz(idx, idy, idz); ...
                     dBy_dy(idx, idy, idz); ...
                     dBy_dz(idx, idy, idz); ...
                     BxLocal(idx, idy, idz); ...
                     ByLocal(idx, idy, idz); ...
                     BzLocal(idx, idy, idz)];

                % Construct Moment matrix
                mx = Moment(1);
                my = Moment(2);
                mz = Moment(3);
                MagMoments = [mx, my, mz,   0,  0,   0,   0,  0
                               0, mx, 0,   my, mz,   0,   0,  0
                             -mz,  0, mx, -mz, my,   0,   0,  0
                               0,  0, 0,    0,  0,   0, -mz,  my
                               0,  0, 0,    0,  0,  mz,   0, -mx
                               0,  0, 0,    0,  0, -my,  mx,  0];

                % Apply the maths to get [Fx;Fy;Fz;Taux;Tauy;Tauz] 
                obj.MagForceTorque(:,i) = MagMoments*U;
            end


        end

        % Function to find the closest grid point
        function [idx, idy, idz] = findClosestGridPoint(obj, Pos)
            [~, idx] = min(abs(obj.x(1,:,1) - Pos(1)));
            [~, idy] = min(abs(obj.y(:,1,1) - Pos(2)));
            [~, idz] = min(abs(obj.z(1,1,:) - Pos(3)));
        end 

        % Function to evalutate gravity affect on each link
        function obj = EvaluateGravity(obj)

            %Get the tentacle Link z positions
            LinkPos = obj.Tentacle.getLinks();
            LinkPosZ = LinkPos(3,:);

            % Get Link Mass
            LinkMass = obj.Tentacle.getLinkMass();

            % Work out height reference system, by finding the lowest
            % possible point in the system
            LinkLength = sqrt((LinkPos(1,1)-LinkPos(1,2))^2+(LinkPos(2,1)-LinkPos(2,2))^2+(LinkPos(3,1)-LinkPos(3,2))^2);
            NumLinks = size(LinkPos,2);
            TotalLength = (1/2)*LinkLength + (NumLinks-1)*LinkLength;
            zBaseline = 0-TotalLength;

            %Evaluate Height from baseline
            H = LinkPosZ - zBaseline;

            % Evaluate Force Due to Gravity
            obj.Fg = H*obj.Gravity*LinkMass;

            %set any tiny elements to 0
            obj.Fg(abs(obj.Fg) < (1.0e-10)) = 0;
        end

        %function to evaluate the effects of joint stiffness
        function obj = EvaluateStiffnes(obj)

            % Get the joint angles
            angles = obj.Tentacle.getJointAngles;

            % Obtain only the alpha values, which is what we are actually
            % interested in
            angles = angles(:,2);

            % Obtain the stiffness value
            stiffness = obj.Tentacle.getStiffness();

            % Apply restoring force due to the stiffness
            % #TODO revisit this method
            restoringTorques = -(stiffness.*angles);

            % Remove the final value since we will ba applying the moment
            % as a force at the COM of each link (fence post, fence panel
            % problem). Ie, the torque at joint 1 (origin) will be applied
            % at the Link center of link 1.
            %   0----x----0----x----0----
            restoringTorques(end) = [];

            %Get the tentacle Link positions and therefore lengths
            LinkPos = obj.Tentacle.getLinks();
            LinkLength = sqrt((LinkPos(1,1)-LinkPos(1,2))^2+(LinkPos(2,1)-LinkPos(2,2))^2+(LinkPos(3,1)-LinkPos(3,2))^2);

            % Evaluate the force to apply to the link COM
            restoringForces = restoringTorques.*LinkLength;


        end

    end
end

