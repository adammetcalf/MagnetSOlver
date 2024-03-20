classdef Tentacle
    % This is a tentacle composed of joints 
    % #TODO update this description
    % after moment intergration etc
    
    %% Private Properties
    properties (Access = private)
        Joints Joint; % An array of joint objects
        HGMs double; % A matrix of the homogenous transfomration matrices
        Jacobean double; % the Jacobean matrix for this tentacle.
        Links double; %COM of each Link
        MagneticMoments double; % A Matrix of magnetic moments Associated with each Link.
        MomentStrength double; % The magnetic moment magnitude for each Link.
    end

    %% Public Methods
    methods

        %% Constructor
        % Constructor
        function obj = Tentacle(LinkLength,Angles,Magnetisation)
            %TENTACLE Construct an instance of Tentacle
            
            % Create n Joints, where n is the number of rows in the
            % provided 'Angles'
            for i=1:size(Angles,1)

                if i == 1
                    % Joint1 is at origin, has no length
                    joint = Joint(Angles(i,1),Angles(i,2),0,0);
                    obj.Joints(i) = joint;
                
                else
                   % Lengths now involved
                   joint = Joint(Angles(i,1),Angles(i,2),0,LinkLength);
                   obj.Joints(i) = joint;
                end
            end

            % Obtain Moment Magnitude of each Link.
            obj.MomentStrength = Magnetisation/(size(Angles,1)-1);

            % Evaluate Homogeneous Transformation Matrices
            obj = EvaluateHGM(obj);

            % Evaluate Jacobean Matrix
            obj = EvaluateJacobean(obj);

            % Evaluate Magnetic Moment of each joint
            obj = EvaluateMagMoments(obj);

        end
        % End Constructor

        %% Update Angles Function
        % Update angles
        function obj = UpdateAngles(obj, Angles)

            if size(Angles,1) ~= length(obj.Joints)
                disp("Incorrect number of angles")
            else
                %Update the angles
                for i=1:size(Angles,1)
                    joint = obj.Joints(i);  %Extract joint
                    obj.Joints(i) = joint.UpdateAngles(Angles(i,:)); %Inject new angle and update joints array
                end
                %End For Loop
            end
            %end If/Else statement

            % Evaluate Homogeneous Transformation Matrices
            obj = EvaluateHGM(obj);
    
            % Evaluate Jacobean Matrix
            obj = EvaluateJacobean(obj);

            % Evaluate Magnetic Moment of each joint
            obj = EvaluateMagMoments(obj);
    
            %End Update Angles
            end

        %% Accessors
        function HGMs = getHGMs(obj)
            % Accessor to retrieve the private property HGMs
            HGMs = obj.HGMs;
        end

        function Jacobean = getJacobean(obj)
            % Accessor to retrieve the private property HGMs
            Jacobean = obj.Jacobean;
        end

        function Links = getLinks(obj)
            % Accessor to retrieve the private property HGMs
            Links = obj.Links;
        end

        function MagneticMoments = getMagneticMoments(obj)
            % Accessor to retrieve the magnetic moment vectors
            MagneticMoments = obj.MagneticMoments;
        end
    end

    %% Private methods
    methods (Access = private)
        
        % Evaluate the Homogeneous Transformation Matrices
        function obj = EvaluateHGM(obj)
            
            % Define origin
            Origin = [1 0 0 0
                      0 1 0 0
                      0 0 1 0
                      0 0 0 1];

            for i=1:length(obj.Joints)

                % Get Joint
                J = obj.Joints(i);
                
                % Extract Frame from Joint
                Frame = getFrame(J);

               
                if i == 1
                % initialise using transform from origin to frame 1

                    % Evaluate homongeneous transformation matrix
                    HGM = Origin*Frame;
                    obj.HGMs(:,:,i) = HGM;

                
                else
                % Transform from previous HGM to next HGM
                    
                    % Evaluate homongeneous transformation matrix
                    HGM = HGM*Frame;
                    obj.HGMs(:,:,i) = HGM;

                end
            %end for loop
            end
        
        % End EvaluateHGM Function
        end

        % Function to evaluate the Jacobean Matrix
        function obj = EvaluateJacobean(obj)
        % The inputs to this function are the homogenous transformation matrices,  the output is the Jacobian.
              
            % Get the number of columns
            NumCols = length(obj.Joints);

            IDvec = [0;0;1];

            %Init Jacobean
            obj.Jacobean = zeros(6,NumCols-1);

            for i = 1:(length(obj.Joints)-1)
                
                    obj.Jacobean(1:3,i) = cross((obj.HGMs(1:3,1:3,i)*IDvec),(obj.HGMs(1:3,4,NumCols)-obj.HGMs(1:3,4,i)));
                    obj.Jacobean(4:6,i) = obj.HGMs(1:3,1:3,i)*IDvec;

            end

            %set any tiny elements to 0
            obj.Jacobean(abs(obj.Jacobean) < (1.0e-10)) = 0;

        end

        % Function to evaluate the magentic moment (magnitude and
        % direction and location) for each Link
        function obj = EvaluateMagMoments(obj)
            
            %Get COMs of each link for applying the moment
            LinkCOMFrames = GetCOM(obj);

            %for each linkCOMframe, work out moment
            for i = 1:size(LinkCOMFrames,3)
                
                % Obtain HGM
                frame = LinkCOMFrames(:,:,i);

                % Populate links matrix
                obj.Links(:,i) = frame(1:3,4);

                % Obtain rotation
                Rotation = frame(1:3, 1:3);

                % Define magnetic moment in local frame orientation
                momentLocal = [0;0;obj.MomentStrength];

                % Transform the magnetic moment to the global frame
                obj.MagneticMoments(:,i) = Rotation * momentLocal;

            end


        end

        % Function to evaluate the COM of each Link
        function LinkCOMFrames = GetCOM(obj)

            % extract start and end HGMs
            for i=1:(length(obj.Joints)-1)

                % Extract Positions (translation components)
                P1 = obj.HGMs(1:3,4,i);
                P2 = obj.HGMs(1:3,4,i+1);

                % Calculate Midpoint
                midPoint = (P1 + P2) / 2;

                % Copy start frame (for correct orientation)
                LinkCOMFrames(:,:,i) = obj.HGMs(:,:,i);

                %Update Postion
                LinkCOMFrames(1:3,4,i) = midPoint;

            end

        end

    end

%% End Class Definition    
end

