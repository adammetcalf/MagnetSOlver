function fig = simpleApp
% SIMPLEAPP Interactively explore plotting functions
%   Choose the function used to plot the sample data to see the
%   differences between surface plots, mesh plots, and waterfall plots

% Create the global variable to store the result
global result;

% Create figure window
fig = uifigure;
fig.Name = "Design Tentacle";

% Manage app layout
gl = uigridlayout(fig,[3 4]);  % 3 rows, 2 columns
gl.RowHeight = {30, '1x', 30}; % Added height for the new button row
gl.ColumnWidth = {'fit', '1x'};

% Create UI components
NumLinks = uieditfield(gl, "numeric");      % To hold the number of links
NumLinksTxt = uilabel(gl);                  % To hold the label for number of links
LinkLength = uieditfield(gl, "numeric");      % To hold the number of links
LinkLengthTxt = uilabel(gl);                  % To hold the label for number of links

NumLinks.Value = 1; % Always at least 1 link
NumLinks.Limits = [1 15]; % Limits from 1 to 15 links
NumLinks.RoundFractionalValues = 'on'; % Always integer
NumLinksTxt.Text = "Number Of Links: ";

LinkLength.Value = 0.01;
LinkLength.Limits = [0.001 0.1];
LinkLengthTxt.Text = "Link Length: ";

ax = uiaxes(gl);
btnClose = uibutton(gl, 'push');  % Add a push button

% Lay out UI components
NumLinksTxt.Layout.Row = 1;
NumLinksTxt.Layout.Column = 1;
NumLinks.Layout.Row = 1;
NumLinks.Layout.Column = 2;
LinkLengthTxt.Layout.Row = 1;
LinkLengthTxt.Layout.Column = 3;
LinkLength.Layout.Row = 1;
LinkLength.Layout.Column = 4;
ax.Layout.Row = 2;
ax.Layout.Column = [1 4];
btnClose.Layout.Row = 3;
btnClose.Layout.Column = [1 2];

% Configure UI component appearance
surf(ax, peaks);
btnClose.Text = 'Close';

% Assign callback function to drop-down
%dd.ValueChangedFcn = {@changePlotType, ax};

% Assign callback to button
btnClose.ButtonPushedFcn = {@closeApp, fig};

    function closeApp(src, event, fig)
        % Store result in the global variable
        result = 'finished';
        delete(fig); % Close the figure and delete it
    end
end

% Program app behavior
function changePlotType(src, event, ax)
    type = event.Value;
    switch type
        case "Surf"
            surf(ax, peaks);
        case "Mesh"
            mesh(ax, peaks);
        case "Waterfall"
            waterfall(ax, peaks);
    end
end
