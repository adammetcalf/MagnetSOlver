function fig = GAapp

% Create figure window
fig = uifigure;
fig.Name = "Current Best Tentacle";

% Manage app Layout
gl = uigridlayout(fig,[2 6]);  % 4 rows, 6 columns
gl.RowHeight = {30,'4x'};  % 
gl.ColumnWidth = {'1x', '1x', '1x', '1x', '1x' };

% Create UI components
EpochTxtLabel = uilabel(gl, 'Text','Epoch: ');                  % To hold the label for number of links
EpochLabel = uilabel(gl);                                       % To hold the Epoch Value
EvoTxtLabel = uilabel(gl, 'Text','Evolution: ');                % To hold the label for number of links
EvoLabel = uilabel(gl);                                         % To hold the Evolution Value
FitnessTxtLabel = uilabel(gl,'Text','Fitness: ');               % To hold the label for fitness
FitnessLabel = uilabel(gl);                                     % To hold the fitness value
ax = uiaxes(gl);                                                % To hold the graph


% Layout UI components
EpochTxtLabel.Layout.Row = 1;
EpochTxtLabel.Layout.Column = 1;
EpochLabel.Layout.Row = 1;
EpochLabel.Layout.Column = 2;
EvoTxtLabel.Layout.Row = 1;
EvoTxtLabel.Layout.Column = 3;
EvoLabel.Layout.Row = 1;
EvoLabel.Layout.Column = 4;
FitnessTxtLabel.Layout.Row = 1;
FitnessTxtLabel.Layout.Column = 5;
FitnessLabel.Layout.Row = 1;
FitnessLabel.Layout.Column = 6;

ax.Layout.Row = 2;
ax.Layout.Column = [1 6];


end