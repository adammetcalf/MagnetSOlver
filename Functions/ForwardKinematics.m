function Positions = ForwardKinematics(DHFrames,Angles)

    % This function performs the forward kinematics. Input 1 is a 3D matrix,
    % where each subsequent DH frame is found at index i from the form
    % (:,:,i). Input 2 is an array of angles, where the row number
    % corresponds the the frame number: [theta1,alpha1;theta2,alpha2.....]


    %get number of frames
    NumFrames = size(DHFrames,3)

    for i = 1:NumFrames
        Frame = DHFrames(:,:,i)

    end


end