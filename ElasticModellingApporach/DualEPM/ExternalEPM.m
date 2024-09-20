classdef ExternalEPM
    properties (Access = private)
        EPMLocal = [0;0;970.1];     % Local magnetic moment (always in z direction)
        Frame double;               % 4x4 transformation matrix of the EPM
        EPMMagneticMoment double;   % Magnetic moment in global coordinates
    end

    methods
        %% Constructor
        function obj = ExternalEPM(location, orientation)
            obj.Frame = eye(4);
            obj = obj.updatePose(location, orientation);
        end

        %% Update Pose Method
        function obj = updatePose(obj, location, orientation)
            obj.Frame(1:3,1:3) = orientation;
            obj.Frame(1:3,4) = location;
            obj.Frame(4,:) = [0,0,0,1];
            obj = obj.EvaluateMagMoment();
        end

        %% Accessor Methods
        function [Position, Moment] = getMag(obj)
            Position = obj.Frame(1:3,4);
            Moment = obj.EPMMagneticMoment;
        end
    end

    methods (Access = private)
        function obj = EvaluateMagMoment(obj)
            Rotation = obj.Frame(1:3,1:3);
            obj.EPMMagneticMoment = Rotation * obj.EPMLocal;
        end
    end
end