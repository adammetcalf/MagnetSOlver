classdef Tentacle
    % This is a tentacle composed of joints 
    % #TODO update this description
    % after moment intergreation etc
    
    %% Private Properties
    properties
        Joints Joint; % An array of joint objects
        HGMs double; % An matrix of the homogenous transfomration matrices
    end

    %% Public Methods
    methods

        % Constructor
        function obj = Tentacle(LinkLength,Angles)
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

            % Evaluate Homogeneous Transformation Matrices
            obj = EvaluateHGM(obj);
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

    end

%% End Class Definition    
end

