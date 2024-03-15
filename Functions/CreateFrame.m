function DHFrame = CreateFrame(theta,alpha,d,a)

    % Create a standard Denavit Hartenberg Frame using theta, alpha, d
    % length and a length inputs.Angles should be input as degrees, and
    % lengths in m

    % convert to radians
    theta = deg2rad(theta);
    alpha = deg2rad(alpha);


    % create frame
    DHFrame = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta)
             sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta)
             0 sin(alpha) cos(alpha) d
             0 0 0 1];

    %set any tiny elements to 0
    DHFrame(DHFrame < (1.0e-10)) = 0;

end