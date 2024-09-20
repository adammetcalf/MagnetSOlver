classdef RobotWithMagnet
    %RobotWithMagnet Class to represent a robot with an active magnet at the end effector

    properties (Access = private)
        robotModel      % Rigid body tree of the robot
        ikSolver        % Inverse kinematics solver
        externalEPM ExternalEPM;    % ExternalEPM object representing the magnet
        currentConfig   % Current joint configuration
        jointLimits     % Joint limits (in radians)
    end

    %% Public Methods
    methods (Access = public)

        %% Constructor
        function obj = RobotWithMagnet(urdfFile)
            
            % load robot
            obj.robotModel = importrobot(urdfFile, 'DataFormat', 'row');

            % Initialize inverse kinematics solver
            obj.ikSolver = generalizedInverseKinematics('RigidBodyTree', obj.robotModel, ...
                'ConstraintInputs', {'position', 'aiming', 'joint'});

            % Initialize ExternalEPM
            obj.externalEPM = ExternalEPM([0;0;0], eye(3));

            % Set joint limits
            obj.jointLimits =  [-170,170;
               -90,120; % Elbow Up
               -170,170;
               -120,120;
               -170,170;
               -120,120;
               -175,175];

            % Set initial configuration
            obj.currentConfig = obj.robotModel.homeConfiguration;

            % Update magnet pose based on initial configuration
            obj = obj.updateMagnetPose();
        end

        %% Update Configuration Method
        function obj = updateConfiguration(obj, jointAngles)

            % Update the robot's joint configuration
            obj.currentConfig = jointAngles;

            % Update the magnet's pose
            obj = obj.updateMagnetPose();
        end

        %% Compute Inverse Kinematics Method
        function [jointAngles, success] = computeIK(obj, targetPosition, targetPoint)
            
            % Compute inverse kinematics to reach a target position and aiming point
            % Define position constraint
            posConstraint = constraintPositionTarget('magnet_center_link');
            posConstraint.TargetPosition = targetPosition;

            % Define aiming constraint
            aimConstraint = constraintAiming('magnet_center_link');
            aimConstraint.TargetPoint = targetPoint;

            % Define joint bounds constraint
            jointConstraint = constraintJointBounds(obj.robotModel);
            jointConstraint.Bounds = obj.jointLimits;

            % Perform inverse kinematics
            [jointAngles, solutionInfo] = obj.ikSolver(obj.currentConfig, posConstraint, aimConstraint, jointConstraint);
            success = solutionInfo.ExitFlag == 1;

            if success
                % Update configuration and magnet pose
                obj = obj.updateConfiguration(jointAngles);
            else
                warning('Inverse kinematics did not converge.');
            end
        end

        %% Show Robot Method
        function showRobot(obj)

            % Visualize the robot in its current configuration
            show(obj.robotModel, obj.currentConfig, 'Frames', 'off', 'PreservePlot', false);
        end

        %% Get Magnetic Field Method
        function [position, moment] = getMagneticField(obj)

            % Get the magnet's position and magnetic moment
            [position, moment] = obj.externalEPM.getMag();
        end

        %% Accessors
        function currentConfig = getCurrentConfig(obj)

            currentConfig = obj.currentConfig;
        end

        function jointLimits = getJointLimits(obj)
            
            jointLimits = obj.jointLimits;
        end
        
    end

    %% Private Methods
    methods (Access = private)
        %% Update Magnet Pose Method
        function obj = updateMagnetPose(obj)
            % Update the magnet's position and orientation based on current configuration
            tform = getTransform(obj.robotModel, obj.currentConfig, 'magnet_center_link', 'base_link');
            position = tform(1:3,4);
            rotation = tform(1:3,1:3);
            % Update ExternalEPM object
            obj.externalEPM = obj.externalEPM.updatePose(position, rotation);
        end
    end
end
