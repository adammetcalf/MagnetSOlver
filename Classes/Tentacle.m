classdef Tentacle
    % This is a tentacle composed of joints and magnetic moments. The
    % tentacle positions (joints and links) is evaluated using homogeneous
    % transfomration matrices and also contains the physical
    % parameters/material properties
    
    %% Private Properties
    properties (Access = private)
        Joints; % An array of joint objects
        HGMs double; % A matrix of the homogenous transfomration matrices
        LinkHGMs double % A matrix of the homogeneous transformation matrices for each link.
        Jacobean double; % the Jacobean matrix for this tentacle.
        Links double; %COM of each Link
        MagneticMoments double; % A Matrix of magnetic moments Associated with each Link.
        MagDirections double; % A Matrix of local magnetic moment directions associated with each link
        MomentStrength double; % The magnetic moment magnitude for each Link.
        JointStiffness double; % Holds the joint stiffness value (single value, not matrix)
        LinkMass double; % holds the mass of each link
        Angles double; %Holds the joint angles
        LinkLength double; %Holds the link length
        Dimensions = 2; % Holds the dimension in which this tentacle may be analysed
    end

    %% Public Methods
    methods

        %% Constructor
        % Constructor
        function obj = Tentacle(LinkLength,Angles,Magnetisation,MagDirections)
            %TENTACLE Construct an instance of Tentacle

            obj.Joints = cell(1);

            obj.LinkLength = LinkLength;

            obj.Angles = Angles;

            obj.MagDirections = MagDirections;
            
            if size(Angles,2) ==2  % 2 dimensional analysis

                obj.Dimensions = 2;

                % Create n Joints, where n is the number of rows in the
                % provided 'Angles'
                for i=1:size(Angles,1)
    
                    % Note: A joint is defined by a denavit hartenberg frame.
                    %joint = Joint(theta,alpha,a,d)
    
                    if i == 1
                        % Joint1 is at origin, has no length
                        joint = Joint(Angles(i,1),Angles(i,2),0,0);
                        obj.Joints{i} = joint;
                    
                    else
                       % Lengths now involved
                       joint = Joint(Angles(i,1),Angles(i,2),0,LinkLength);
                       obj.Joints{i} = joint;
                    end
                end

            else % 3 dimensional analysis

                obj.Dimensions = 3;

                % Create n Joints, where n is the number of rows in the
                % provided 'Angles'
                for i=1:size(Angles,1)
    
                    % Note: A joint is defined by a denavit hartenberg frame.
                    %joint = Joint(theta,alpha,a,d)
    
                    if i == 1
                        % Joint1 is at origin, has no length
                        joint = Joint3D(Angles(i,1),Angles(i,2),Angles(i,3),0,0);
                        obj.Joints{i} = joint;
                    
                    else
                       % Lengths now involved
                       joint = Joint3D(Angles(i,1),Angles(i,2),Angles(i,3),0,LinkLength);
                       obj.Joints{i} = joint;
                    end
                end      

            end

            % Link Radii (m)
            radius = 1e-03;
            
            % EcoFlex 0030 Density (kg m^-3)
            rho = 1070; 
            
            % Youngs Modulus Ecoflex 0030 (125 kPa)
            E = 125000;

            % Evaluate Material properties
            obj = MaterialProperties(obj, radius,rho,LinkLength,E);

            % Obtain Moment Magnitude of each Link.
            obj.MomentStrength = Magnetisation/(size(Angles,1)-1);

            % Evaluate Homogeneous Transformation Matrices of joints
            obj = EvaluateHGM(obj);

            % Evaluate of Homogeneous Transformation Matrices links
            obj = GetCOM(obj);

            % Evaluate Jacobean Matrix
            obj = EvaluateJacobean(obj);

            % Evaluate Magnetic Moment of each joint
            obj = EvaluateMagMoments(obj);

        end
        % End Constructor

        %% Update Angles Function
        % Update angles
        function obj = UpdateAngles(obj, Angles)

            obj.Angles = Angles;

            if size(Angles,1) ~= length(obj.Joints)
                disp("Incorrect number of angles")
            else
                %Update the angles
                for i=1:size(Angles,1)
                    joint = obj.Joints{i};  %Extract joint
                    obj.Joints{i} = joint.UpdateAngles(Angles(i,:)); %Inject new angle and update joints array
                end
                %End For Loop
            end
            %end If/Else statement

            % Evaluate Homogeneous Transformation Matrices
            obj = EvaluateHGM(obj);
     
            % Evaluate of Homogeneous Transformation Matrices links
            obj = GetCOM(obj);
    
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
            % Accessor to retrieve the private property link positions
            Links = obj.Links;
        end

        function MagneticMoments = getMagneticMoments(obj)
            % Accessor to retrieve the magnetic moment vectors
            MagneticMoments = obj.MagneticMoments;
        end

        function JointStiffness = getStiffness(obj)
            % Accessor to retrieve the Joint stiffness
            JointStiffness = obj.JointStiffness;
        end

        function LinkMass = getLinkMass(obj)
            % Accessor to retrieve the Joint stiffness
            LinkMass = obj.LinkMass;
        end

        %get the joint angles
        function Angles = getJointAngles(obj)
            Angles = obj.Angles;
        end

        % get the link HGMS
        function LinkHGMs = getLinkHGMs(obj)
            LinkHGMs = obj.LinkHGMs;
        end

        % get the link Length
        function LinkLength = getLinkLength(obj)
            LinkLength = obj.LinkLength;
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
                J = obj.Joints{i};
                
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
            
            %for each linkCOMframe, work out moment
            for i = 1:size(obj.LinkHGMs,3)
                
                % Obtain HGM
                frame = obj.LinkHGMs(:,:,i);

                % Populate links matrix
                obj.Links(:,i) = frame(1:3,4);

                % Obtain rotation
                Rotation = frame(1:3, 1:3);

                % Obtain local magnetic moment
                momentLocal = obj.MagDirections(:,i)*obj.MomentStrength;

                % Transform the magnetic moment to the global frame
                obj.MagneticMoments(:,i) = Rotation * momentLocal;

            end


        end

        % Function to evaluate the COM of each Link
        function obj = GetCOM(obj)

            % extract start and end HGMs
            for i=1:(length(obj.Joints)-1)

                % Extract Positions (translation components)
                P1 = obj.HGMs(1:3,4,i);
                P2 = obj.HGMs(1:3,4,i+1);

                % Calculate Midpoint
                midPoint = (P1 + P2) / 2;

                % Copy start frame (for correct orientation)
                obj.LinkHGMs(:,:,i) = obj.HGMs(:,:,i);

                %Update Postion
                obj.LinkHGMs(1:3,4,i) = midPoint;

            end

        end

        % obtain link mass and joint stiffness
        function obj = MaterialProperties(obj, radius,rho,Length,E)

            % Link Volumes (m^3)
            Vol = Length*pi*(radius^2);       
            
            % Link Masses 
            obj.LinkMass = rho*Vol;
            
            %%%%% Evalaute inertia and stiffness using link length
            
            % Define Inertia (kg m^-2) 
            %%% TO DO: Inertia evaluated about correct position?
            %%% SOLID CYCLNDER ABOUT ENDPOINT DIAMETER
            %%% I = 1/4*m*R^2 + 1/3*m*L^2
            
            %I = (1/4)*obj.LinkMass*radius^2 + (1/3)*obj.LinkMass*Length;

            %obj.JointStiffness = (E*I)/Length;
            %Note, the above is for longitudinal stiffness (ie, resistance
            %to stretch)

            % Ending stiffness calculation
            I  = (1/4)*pi*(radius^4);
            obj.JointStiffness = (E*I);

        end

    end

%% End Class Definition    
end

