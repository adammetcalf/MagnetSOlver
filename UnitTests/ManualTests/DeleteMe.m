close all;
clear;
clc;

% Initialise tentacle and angles
angles = [90,90;0,10;0,20;0,30;0,40;0,50];
Magdirection = [0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
TentacleMagnetisation = 12;
LinkLength = 0.01;

%Get initial angle
initalAlpha = angles(1,2);

% Add 10 to initial alpha just so we have an obvious difference
initalAlpha = initalAlpha+15;

% spin up a cheeky tentacle
Tentacle = Tentacle(LinkLength,angles,TentacleMagnetisation,Magdirection);

%clear certain variables just for workspace navigation ease
clear Magdirection TentacleMagnetisation LinkLength

% Obtain only the alpha values, which is what we are actually
% interested in
angles = angles(:,2);

% Joint 1 (at the origin) must be offset from the initial
% value, since the initial value is used to orientate the
% tentacle at its 'rest' position
angles(1) = initalAlpha - angles(1); % #TODO - account for 180 to -180 overlap?

% Obtain the stiffness value
stiffness = Tentacle.getStiffness();

% Apply restoring force due to the stiffness
% #TODO revisit this method
restoringTorques = -(stiffness.*angles);

% Remove the final value since we will ba applying the moment
% as a force at the COM of each link (fence post, fence panel
% problem). Ie, the torque at joint 1 (origin) will be applied
% at the Link center of link 1.
%   0----x----0----x----0----
restoringTorques(end) = [];

%Get the tentacle Link positions and therefore lengths
LinkPos = Tentacle.getLinks();
LinkLength = sqrt((LinkPos(1,1)-LinkPos(1,2))^2+(LinkPos(2,1)-LinkPos(2,2))^2+(LinkPos(3,1)-LinkPos(3,2))^2);

% Evaluate the force to apply to the link COM (which is half of
% a linklength away from the joint)
restoringForces = restoringTorques.*(LinkLength/2);

% Get Joint Frames
Jointframes = Tentacle.getHGMs();

% Get Link Frames
Linkframes = Tentacle.getLinkHGMs();

% Assuming restoringForces are correctly calculated and represent forces in the local z-direction of each link

% Initialize global force vectors
globalForces = zeros(size(Linkframes,3), 3); % Assuming Linkframes is a 4x4xN 3D matrix, where N is the number of links

for i = 1:size(Linkframes,3)
    
    localForceDirection = [0; 1; 0]; % Update this based on actual force direction
    localForceMagnitude = restoringForces(i); % The calculated force magnitude
    localForce = localForceDirection * localForceMagnitude;
    
    % Homogeneous representation of local force
    localForceHom = [localForce; 0]; % Append 0 for homogeneous representation
    
    % Transform the force vector to global coordinates using the link's HTM
    globalForceHom = Linkframes(:,:,i) * localForceHom;
    
    % Extract global force vector (discard the homogeneous coordinate)
    globalForce(:,i) = globalForceHom(1:3);
    
    % Apply the global force vector as needed in your model
end



