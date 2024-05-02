classdef Joint < JointInterface
    % This is a joint, which is represented by a single standard Denavit
    % Hartenberg frame
    
    %% Private Properties
    properties (Access = private) 
        a double; %DH Parameter (in m)
        d double; %DH Parameter (in m)
        theta double; %DH Parameter (in radians)
        alpha double; %DH Parameter (in radians)
        alpha2 double; %DH Parameter (in radians) -- NOT USED IN THIS CHILD< BUT NECESSARY FOR INHERITENCE
        Frame double = zeros(4,4); % DH Frame
    end
    
    %% Public Methods
    methods
        
        %Constructor
        function obj = Joint(theta,alpha,a,d)
            %JOINT Construct an instance of Joint

            % Inject DH parameters
            obj.a = a;
            obj.d = d;
            obj.theta = deg2rad(theta);
            obj.alpha = deg2rad(alpha);

            % Evaluate Frame
            obj = evaluateFrame(obj);
        end

        function obj = UpdateAngles(obj,Angles)
            % Update the angles and re-evaluate the frame

            obj.theta = deg2rad(Angles(1));
            obj.alpha = deg2rad(Angles(2));

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

            obj.Frame = [cos(obj.theta),-sin(obj.theta)*cos(obj.alpha),sin(obj.theta)*sin(obj.alpha),obj.a*cos(obj.theta);...
            sin(obj.theta),cos(obj.theta)*cos(obj.alpha),-cos(obj.theta)*sin(obj.alpha),obj.a*sin(obj.theta);...
            0,sin(obj.alpha),cos(obj.alpha),obj.d;...
            0,0,0,1];

            %set any tiny elements to 0
            obj.Frame(abs(obj.Frame) < (1.0e-10)) = 0;
        end
    
    end

%% End Class definition
end

