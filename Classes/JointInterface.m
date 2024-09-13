classdef JointInterface
    %JOINTINTERFACE parent of Joint and Joint3D
    

    
    %% Abstract Properties
    properties (Abstract, Access = protected)
        a double;      % DH Parameter (in m)
        d double;      % DH Parameter (in m)
        theta double;  % DH Parameter (in radians) - corresponding to local z deflection
        alpha double;  % DH Parameter (in radians) - corresponding to local x deflection
        alpha2 double; % DH Parameter (in radians) - corresponding to local y deflection
        Frame double;  % DH Frame
    end

    
    
    %% Public Methods
    methods
        obj = UpdateAngles(obj, Angles) % Method to update joint angles
        Frame = getFrame(obj)           % Accessor to retrieve the DH frame
    end

    %% Privste methods
    methods (Access = private)
        obj = evaluateFrame(obj)        % Evaluates the frame assoicated with the joint
    end
    
end
   


