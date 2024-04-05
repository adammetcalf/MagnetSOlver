classdef ExternalEPM
    %EXTERNALEPM This class introduces an external magnet, which can eb
    %added to the field effects of the world class
    
    
    properties (Access = private)
        EPMLocal = [0;0;970.1]      % Define the local orientation of the moment (always in z direction)
        Frame double;               % A 4x4 matrix represtning the position and orientation of the EPM
        EPMMagneticMoment double;   % Magnetic moment in global.
    end
    
    %% Public Methods
    methods (Access = public)

        %% Constructor
        function obj = ExternalEPM(location,orientation)
            % Uses the position and orientation of the EPM to construct a
            % frame
            obj.Frame(1:3,1:3) = orientation;
            obj.Frame(1:3,4) = location;
            obj.Frame(4,:) = [0,0,0,1];

            obj =  EvaluateMagMoment(obj);
        end
        
        %% Accessors
        function [Position, Moment] = getMag(obj)

            %Get the position vector from the frame
            Position = obj.Frame(1:3,4);

            Moment = obj.EPMMagneticMoment;
        end

        
    end

    %% Private Methods
    methods (Access = private)

        function obj = EvaluateMagMoment(obj)

            % Get the rotation of the the EPM frame
            Rotation = obj.Frame(1:3,1:3);

            % Transform the magnetic moment to the global frame
            obj.EPMMagneticMoment = Rotation * obj.EPMLocal;

        end



    end
end

