function fig = simpleApp
% SIMPLEAPP 

% Create the global variable to store the result
global tentacleResult;

%% Create UI Components

% Create figure window
fig = uifigure;
fig.Name = "Design Tentacle";

% Manage app layout
gl = uigridlayout(fig,[4 5]);  % 5 rows, 5 columns
gl.RowHeight = {30, 150, '4x', 30};  % Set second row to fixed height for the table
gl.ColumnWidth = {'1x', '1x', '1x', '1x', '1x' };

% Create UI components
NumLinks = uieditfield(gl, "numeric","ValueChangedFcn",@(src,event) NumLinksValueChanged(src,event,fig));                                  % To hold the number of links
NumLinksTxt = uilabel(gl, 'Text','Number Of Links: ');                  % To hold the label for number of links
LinkLength = uieditfield(gl, "numeric","ValueChangedFcn",@(src,event) LengthValueChanged(src,event,fig));                                % To hold the number of links
LinkLengthTxt = uilabel(gl, 'Text','Link Length: ');                    % To hold the label for number of links
btnClose = uibutton(gl, 'push','Text', 'Complete');                     % Add a push button to complete the tentacle creation
btnReset = uibutton(gl, 'push','Text', 'Reset');                        % Add a push button to reset the tentacle creation
ax = uiaxes(gl);                                                        % Axes to hold the tentacle
table = uitable(gl, "CellEditCallback",@(src,event) TableValueChanged(src,event,fig));


% Define the 'number of links' input field
NumLinks.Value = 1;                                                     % Always at least 1 link
NumLinks.HorizontalAlignment = 'center';                                % Centre the field entry
NumLinks.Limits = [1 15];                                               % Limits from 1 to 15 links
NumLinks.RoundFractionalValues = 'on';                                  % Always integer

% Define the 'Link Length' input field
LinkLength.Value = 0.01;                                                % Abritrary Link Length
LinkLength.HorizontalAlignment = 'center';                              % Centre the field entry
LinkLength.Limits = [0.001 0.1];                                        % Limit the lingth lengths

% Define the table
s = uistyle("HorizontalAlignment","center");
table.Data = zeros(4, NumLinks.Value);
table.Data(3:4,:) = 1;
table.ColumnName = arrayfun(@(i) sprintf('Link %d', i), 1:NumLinks.Value, 'UniformOutput', false);
table.RowName = {'X', 'Y', 'Z','NdFeb?'};
for i = 1:NumLinks.Value
    columnEdit(i) = true;
end
table.ColumnEditable = columnEdit;
addStyle(table,s)

% create links 
Angles = zeros(NumLinks.Value+1,2);
Angles(1,:) = [90,90];

tentacleResult = Angles;

% arbitrary magentisation
Magnetisation = 2*NumLinks.Value;

% create tentacle
tentacle = Tentacle(LinkLength.Value,Angles,Magnetisation,table.Data(1:3,:));




%% Lay out UI components
% Row 1
NumLinksTxt.Layout.Row = 1;
NumLinksTxt.Layout.Column = 1;
NumLinks.Layout.Row = 1;
NumLinks.Layout.Column = 2;
LinkLengthTxt.Layout.Row = 1;
LinkLengthTxt.Layout.Column = 4;
LinkLength.Layout.Row = 1;
LinkLength.Layout.Column = 5;

% Row 2
table.Layout.Row = 2;
table.Layout.Column = [1 5];

% Row 3
ax.Layout.Row = 3;
ax.Layout.Column = [1 5];

% Row 4
btnClose.Layout.Row = 4;
btnClose.Layout.Column = [4 5];
btnReset.Layout.Row = 4;
btnReset.Layout.Column = [1 2];

% Configure UI component appearance
plotTentacle(fig);

%% Callbacks

% Assign callback to close button
btnClose.ButtonPushedFcn = {@closeApp, fig};

    function closeApp(src, event, fig)
        delete(fig); % Close the figure and delete it
    end

% Assign callback to reset button
btnReset.ButtonPushedFcn = {@resetApp, fig};

    function resetApp(src, event, fig)
        
        LinkLength.Value = 0.01;
        NumLinks.Value = 1;

        % repopulate table, evaluate everything and create/overwrite tentacle
        table.Data = zeros(4, NumLinks.Value);
        table.Data(3:4,:) = 1;
        table.ColumnName = arrayfun(@(i) sprintf('Link %d', i), 1:NumLinks.Value, 'UniformOutput', false);

        % create links 
        Angles = zeros(NumLinks.Value+1,2);
        Angles(1,:) = [90,90];

        tentacleResult = Angles;
        
        % #TODO sort out the arbitrary magentisation profile
        % arbitrary magentisation
        Magnetisation = 2*NumLinks.Value;
        
        % create tentacle
        tentacle = Tentacle(LinkLength.Value,Angles,Magnetisation,table.Data(1:3,:));

        plotTentacle(fig);

        tentacleResult = tentacle;

    end

% NumLinks value changed callback
    function NumLinksValueChanged(src,event,fig)

        % repopulate table, evaluate everything and create/overwrite tentacle
        table.Data = zeros(4, NumLinks.Value);
        table.Data(3:4,:) = 1;
        table.ColumnName = arrayfun(@(i) sprintf('Link %d', i), 1:NumLinks.Value, 'UniformOutput', false);

        % create links 
        Angles = zeros(NumLinks.Value+1,2);
        Angles(1,:) = [90,90];

        tentacleResult = Angles;
        
        % arbitrary magentisation
        Magnetisation = 2*NumLinks.Value;
        
        % create tentacle
        tentacle = Tentacle(LinkLength.Value,Angles,Magnetisation,table.Data(1:3,:));

        plotTentacle(fig);

        tentacleResult = tentacle;

    end

    function LengthValueChanged(src,event,fig)

        % create links 
        Angles = zeros(NumLinks.Value+1,2);
        Angles(1,:) = [90,90];

        tentacleResult = Angles;
        
        % arbitrary magentisation
        Magnetisation = 2*NumLinks.Value;
        
        % create tentacle
        tentacle = Tentacle(LinkLength.Value,Angles,Magnetisation,table.Data((1:3),:));

        plotTentacle(fig);

        tentacleResult = tentacle;

    end

%function to plot the tentacle
    function plotTentacle(fig)
        % Retrieve the tentacle Joint positions
        HGMs = tentacle.getHGMs();  % Assuming getHGMs returns homogeneous transformation matrices for all joints
        
        % Extract joint positions from the transformation matrices
        jointPos = zeros(3, size(HGMs, 3));
        for j = 1:size(HGMs, 3)
            jointPos(:, j) = HGMs(1:3, 4, j);  
        end
        
        % Retrieve the tentacle Link positions
        LinkPos = tentacle.getLinks();  % Assuming getLinks returns positions for links
    
        % Define the plot extent based on the link length and number of links
        PlotLength = LinkLength.Value * (NumLinks.Value + 2);
    
        % Get magnetic moments
        Moments = tentacle.getMagneticMoments();
    
        % Plotting joint positions
        plot3(ax, jointPos(1, :), jointPos(2, :), jointPos(3, :), ...
              'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'LineWidth', 1);
        hold(ax, 'on');  % Ensure that the plot does not clear existing figures
        
        % Plotting link positions
        plot3(ax, LinkPos(1, :), LinkPos(2, :), LinkPos(3, :), ...
              'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    
        % Plot Moments
        PlotMagenticMoments(fig,LinkPos,Moments);
    
        % Set axes properties
        axis(ax, 'equal');
        axis(ax, [-LinkLength.Value*2 PlotLength -LinkLength.Value*2 LinkLength.Value*2 -LinkLength.Value*2 LinkLength.Value*2]);
        grid(ax, 'on');
        xlabel(ax, 'X');
        ylabel(ax, 'Y');
        zlabel(ax, 'Z');
        hold(ax, 'off');  % Release the hold to allow other plotting commands to function normally
    end

    %Function to plot the magnetic moments
    function PlotMagenticMoments(fig,LinkPos,Moments)

        for j = 1:size(LinkPos, 2)
            Pos = LinkPos(:, j);  % Current link position
            Moment = Moments(:, j);  % Current magnetic moment
    
            % Normalize the moment to have the magnitude of the link length
            NormalizedMoment = (Moment / norm(Moment)) * LinkLength.Value/3;
    
            % Calculate the new starting position to include the z offset
            NewPos = Pos + [0; 0; LinkLength.Value];
    
            % Plot the moment as a quiver, starting from the new position
            q = quiver3(ax, NewPos(1), NewPos(2), NewPos(3), NormalizedMoment(1), NormalizedMoment(2), NormalizedMoment(3), 'k', 'LineWidth', 1);
            q.MaxHeadSize = 4 * q.MaxHeadSize; % Double the size of the arrowhead

            % Plotting link positions
            plot3(ax, NewPos(1, :), NewPos(2, :), NewPos(3, :), ...
            'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
        end

    end

    function TableValueChanged(src,event,fig)

        %Check for material
        CheckNdFeb(fig);

        % create links 
        Angles = zeros(NumLinks.Value+1,2);
        Angles(1,:) = [90,90];
        
        % arbitrary magentisation
        Magnetisation = 2*NumLinks.Value;
        
        % create tentacle
        tentacle = Tentacle(LinkLength.Value,Angles,Magnetisation,table.Data(1:3,:));

        plotTentacle(fig);

        tentacleResult = tentacle;

    end

    function CheckNdFeb(fig)
   % loop though row 4 of the table. If an element == 0, set rows
   % 1,3 to 0

       for j = 1:size(table.Data,1)
           if table.Data(4,j) == 1
              % Do nothing
           else
               table.Data(1:3,j) = 0;
           end

       end

    end


end
