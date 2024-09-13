classdef Joint3D < JointInterface
    % This is a joint, which is represented by a single standard Denavit
    % Hartenberg frame
    
    %% Private Properties
    properties (Access = protected) 
        a double; %DH Parameter (in m)
        d double; %DH Parameter (in m)
        theta double; %DH Parameter (in radians) - corresponding to local z deflection
        alpha double; %DH Parameter (in radians) - corresponding to local x deflection
        alpha2 double;  %DH Parameter (in radians) - corresponding to local y deflection
        Frame double = zeros(4,4); % DH Frame
    end
    
    %% Public Methods
    methods
        
        %Constructor
        function obj = Joint3D(theta,alpha,alpha2,a,d)
            %JOINT Construct an instance of Joint

            % Inject DH parameters
            obj.a = a;
            obj.d = d;
            obj.theta = deg2rad(theta);
            obj.alpha = deg2rad(alpha);
            obj.alpha2 = deg2rad(alpha2);

            % Evaluate Frame
            obj = evaluateFrame(obj);
        end

        function obj = UpdateAngles(obj,Angles)
            % Update the angles and re-evaluate the frame

            obj.theta = deg2rad(Angles(1));
            obj.alpha = deg2rad(Angles(2));
            obj.alpha2 = deg2rad(Angles(3));

            % Evaluate Frame
            obj = evaluateFrame(obj);

        end

        %% Accessors
        function Frame = getFrame(obj)
            
            % Accessor to retrieve the private property DH frame
            Frame = obj.Frame;
        end

    end    
    
    
    %% Private Methods
    methods (Access = private)
        function obj = evaluateFrame(obj)
            %evaluateFrame uses the internal properties to evaluate the DH
            %frame associated with this Joint.

            %Frame 1 (z-x)
            Frame1 = [cos(obj.theta),-sin(obj.theta)*cos(obj.alpha),sin(obj.theta)*sin(obj.alpha),obj.a*cos(obj.theta);...
            sin(obj.theta),cos(obj.theta)*cos(obj.alpha),-cos(obj.theta)*sin(obj.alpha),obj.a*sin(obj.theta);...
            0,sin(obj.alpha),cos(obj.alpha),obj.d;...
            0,0,0,1];

            %Define frame to rotate by 90 degrees about z
            Pos90Z = [0 -1 0 0
                      1 0 0 0
                      0 0 1 0
                      0 0 0 1];

            % Define frame to rotate by -90 degrees about z
            Neg90Z = [0 1 0 0
                      -1 0 0 0
                      0 0 1 0
                      0 0 0 1];

            %Frame 1 (z-y) - rotate frame 1 about z +90, evaluate alpha 2
            %Note that for this transformation, theta, a and d are set to
            %zero because these have already been handled by frame 1

            Frame2 = [1,0,0,0;...
            0,cos(obj.alpha2),-1*sin(obj.alpha2),0;...
            0,sin(obj.alpha2),cos(obj.alpha2),0;...
            0,0,0,1];

            % Perfrom the full set of transformations
            obj.Frame = Frame1*Pos90Z*Frame2*Neg90Z;

            %set any tiny elements to 0
            obj.Frame(abs(obj.Frame) < (1.0e-10)) = 0;
        end
    
    end

%% End Class definition
end

