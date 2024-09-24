function [x_positions, y_positions, moments] = elastica_UI(NumNodes, Length)
    % Parameters
    n = NumNodes;                 % Number of nodes
    L = Length;                   % Total length
    delta_s = L / (n - 1);        % Segment length
    s = linspace(0, L, n)';       % Arc lengths

    % Initial configuration (straight line along x-axis)
    r = zeros(n, 2);
    r(:,1) = linspace(0, L, n);

    % Initialize magnetic moments aligned along x-axis
    moments = repmat([1, 0], n, 1);  % Magnetic moments at each node

    % Create figure and plot
    hFig = figure('Name', 'Interactive Elastica with Magnetic Moments', 'NumberTitle', 'off');
    hAx = axes('Parent', hFig);
    hold(hAx, 'on');
    % Plot magnetic moments using quiver
    hMoments = quiver(hAx, r(:,1), r(:,2), moments(:,1), moments(:,2), 0.5, 'b', 'LineWidth', 1);
    hElastica = plot(hAx, r(:,1), r(:,2), 'k.-', 'LineWidth', 2);
    hPoints = plot(hAx, r(:,1), r(:,2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    axis equal;

    % Set axes limits to +/- (Length + 2 segments)
    xlim([-(L + 2 * delta_s), L + 2 * delta_s]);
    ylim([-(L + 2 * delta_s), L + 2 * delta_s]);

    xlabel('x position (m)');
    ylabel('y position (m)');
    title('Drag points to reshape the elastica (Segment lengths constant)');
    hold off;

    % Create UI controls (buttons) at the bottom row
    hResetButton = uicontrol('Style', 'pushbutton', 'String', 'Reset', ...
        'Units', 'normalized', 'Position', [0.05 0.01 0.1 0.05], 'Callback', @resetButtonCallback);
    hLoadButton = uicontrol('Style', 'pushbutton', 'String', 'Load Jshape', ...
        'Units', 'normalized', 'Position', [0.25 0.01 0.1 0.05], 'Callback', @loadButtonCallback);
    hInvertButton = uicontrol('Style', 'pushbutton', 'String', 'Invert Moments', ...
        'Units', 'normalized', 'Position', [0.45 0.01 0.1 0.05], 'Callback', @invertButtonCallback);
    hStopButton = uicontrol('Style', 'pushbutton', 'String', 'Stop', ...
        'Units', 'normalized', 'Position', [0.8 0.01 0.1 0.05], 'Callback', @stopButtonCallback);

    % Store data for callbacks
    data.r = r;
    data.moments = moments;
    data.hPoints = hPoints;
    data.hElastica = hElastica;
    data.hMoments = hMoments;
    data.delta_s = delta_s;
    data.dragIdx = [];
    data.isRunning = true;
    guidata(hFig, data);

    % Set up draggable points
    set(hPoints, 'ButtonDownFcn', @startDragFcn);
    set(hFig, 'WindowButtonUpFcn', @stopDragFcn);

    % Wait until user presses 'Stop'
    while ishandle(hFig) && data.isRunning
        drawnow;
        data = guidata(hFig);
    end

    % After 'Stop' is pressed, get x_positions and y_positions
    if ishandle(hFig)
        x_positions = data.r(:,1);
        y_positions = data.r(:,2);
        moments = data.moments;
        % Close the figure
        close(hFig);
    else
        x_positions = [];
        y_positions = [];
        moments = [];
    end

    % Callback functions
    function resetButtonCallback(src, ~)
        data = guidata(src);
        % Reset the elastica to horizontal
        data.r = zeros(n, 2);
        data.r(:,1) = linspace(0, L, n);
        % Reset magnetic moments to x-axis
        data.moments = repmat([1, 0], n, 1);
        % Update the plot
        set(data.hPoints, 'XData', data.r(:,1), 'YData', data.r(:,2));
        set(data.hElastica, 'XData', data.r(:,1), 'YData', data.r(:,2));
        set(data.hMoments, 'XData', data.r(:,1), 'YData', data.r(:,2), ...
            'UData', data.moments(:,1), 'VData', data.moments(:,2));
        guidata(src, data);
    end

    function loadButtonCallback(src, ~)
    data = guidata(src);
    % Load Jshape.mat, which should contain a 50x2 matrix
    loaded_data = load('Jshape.mat');
        if isfield(loaded_data, 'Jshape') && size(loaded_data.Jshape, 1) == 50 && size(loaded_data.Jshape, 2) == 2
            % Convert positions from mm to m
            Jshape = loaded_data.Jshape / 1000;
            
            % Smooth the Jshape using a moving average
            windowSize = 5; % Adjust the window size as needed
            Jshape_smooth(:, 1) = smoothdata(Jshape(:, 1), 'movmean', windowSize);
            Jshape_smooth(:, 2) = smoothdata(Jshape(:, 2), 'movmean', windowSize);
            
            % Update the elastica nodes with the smoothed positions
            data.r = Jshape_smooth;
            
            % Update magnetic moments to align with local tangents
            data.moments = compute_moments(data.r);
            
            % Update the plot
            set(data.hPoints, 'XData', data.r(:,1), 'YData', data.r(:,2));
            set(data.hElastica, 'XData', data.r(:,1), 'YData', data.r(:,2));
            set(data.hMoments, 'XData', data.r(:,1), 'YData', data.r(:,2), ...
                'UData', data.moments(:,1), 'VData', data.moments(:,2));
        else
            disp('Jshape.mat must contain a 50x2 matrix.');
        end
        guidata(src, data);
    end

    function invertButtonCallback(src, ~)
        data = guidata(src);
        n = size(data.r, 1);
        % Determine indices for the second half of the tentacle
        idx_start = ceil(n/2) + 1;
        idx_end = n;
        % Invert the moments for the second half
        data.moments(idx_start:idx_end, :) = -data.moments(idx_start:idx_end, :);
        % Update the plot
        set(data.hMoments, 'UData', data.moments(:,1), 'VData', data.moments(:,2));
        guidata(src, data);
    end

    function stopButtonCallback(src, ~)
        data = guidata(src);
        % Set isRunning to false to exit the while loop
        data.isRunning = false;
        guidata(src, data);
    end

    function startDragFcn(src, ~)
        data = guidata(src);
        % Get mouse click position
        pt = get(hAx, 'CurrentPoint');
        cp = pt(1,1:2);
        % Find the closest node to the click
        [~, idx] = min(sum((data.r - cp).^2, 2));
        data.dragIdx = idx;
        guidata(src, data);
        % Start dragging
        set(hFig, 'WindowButtonMotionFcn', @draggingFcn);
    end

    function draggingFcn(src, ~)
        data = guidata(src);
        pt = get(hAx, 'CurrentPoint');
        cp = pt(1,1:2);
        idx = data.dragIdx;
        % Update the position of the dragged node
        data.r(idx, :) = cp;
        % Apply FABRIK to adjust other nodes
        data.r = fabrik(data.r, idx, data.delta_s);
        % Update magnetic moments to align with local tangents
        data.moments = compute_moments(data.r);
        % Update the plot
        set(data.hPoints, 'XData', data.r(:,1), 'YData', data.r(:,2));
        set(data.hElastica, 'XData', data.r(:,1), 'YData', data.r(:,2));
        set(data.hMoments, 'XData', data.r(:,1), 'YData', data.r(:,2), ...
            'UData', data.moments(:,1), 'VData', data.moments(:,2));
        guidata(src, data);
    end

    function stopDragFcn(~, ~)
        set(hFig, 'WindowButtonMotionFcn', '');
    end

    function r_new = fabrik(r, idx, delta_s)
        % Positions of the nodes
        n = size(r, 1);
        % Fixed base node
        base_pos = r(1,:);
        % Target position (position of the dragged node)
        target_pos = r(idx,:);
        % Tolerance for convergence
        tol = 1e-6;
        % Maximum iterations
        max_iter = 100;

        % Initial distance between base and target
        dist = norm(target_pos - base_pos);

        % Check if the target is reachable
        total_length = delta_s * (n - 1);
        if dist > total_length
            % Target is unreachable; stretch towards the target
            for i = 1:n-1
                r(i+1,:) = r(i,:) + (target_pos - r(i,:)) * (delta_s / norm(target_pos - r(i,:)));
            end
        else
            % Target is reachable; proceed with iterations
            b = base_pos;
            diff = inf;
            iter = 0;
            while diff > tol && iter < max_iter
                r_old = r;
                % Backward reaching
                r(idx,:) = target_pos;
                for i = idx-1:-1:1
                    dir = r(i,:) - r(i+1,:);
                    dir = dir / norm(dir);
                    r(i,:) = r(i+1,:) + dir * delta_s;
                end
                for i = idx+1:n
                    dir = r(i,:) - r(i-1,:);
                    dir = dir / norm(dir);
                    r(i,:) = r(i-1,:) + dir * delta_s;
                end
                % Forward reaching
                r(1,:) = b;
                for i = 1:n-1
                    dir = r(i+1,:) - r(i,:);
                    dir = dir / norm(dir);
                    r(i+1,:) = r(i,:) + dir * delta_s;
                end
                % Check for convergence
                diff = max(vecnorm(r - r_old, 2, 2));
                iter = iter + 1;
            end
        end
        r_new = r;
    end

    function moments = compute_moments(r)
        % Compute magnetic moments aligned with local tangents
        n = size(r, 1);
        moments = zeros(n, 2);
        for i = 1:n-1
            dir = r(i+1,:) - r(i,:);
            dir = dir / norm(dir);
            moments(i,:) = dir;
        end
        % For the last node, use the previous segment direction
        moments(n,:) = moments(n-1,:);
    end
end
