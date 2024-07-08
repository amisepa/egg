% CLEAN_SIGNAL Manually select and remove bad segments from a time series.
%
% Usage:
%   [cleanSignal, cleanT] = clean_signal(signal, t)
%
% Inputs:
%   signal - A vector representing the time series data.
%   t      - A vector representing the time points corresponding to the signal.
%
% Outputs:
%   cleanSignal - The cleaned time series with bad segments removed.
%   cleanT      - The adjusted time vector with bad segments removed and made continuous.
%
% Example:
%   % Generate a sample time series
%   t = linspace(0, 10, 1000);
%   signal = sin(2*pi*0.5*t) + 0.5*randn(size(t));
%
%   % Clean the signal by removing bad segments
%   [cleanSignal, cleanT] = clean_signal(signal, t);
%
% Copyright (C), Cedric Cannard, June 2024

function [cleanSignal, cleanT] = clean_signal(signal, t)

    % Initialize figure and plot
    f = figure('Name', 'Manual Signal Cleaning', 'Color', 'white', ...
               'CloseRequestFcn', @closeFigure, 'NumberTitle', 'off', ...
               'MenuBar', 'none', 'ToolBar', 'none');
    ax = axes('Parent', f);
    plot(ax, t, signal, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
    title(ax, 'Select bad segments to remove');
    xlabel(ax, 'Time (minutes)', 'FontWeight', 'bold');
    ylabel(ax, 'Amplitude', 'FontWeight', 'bold');
    
    % Enable interactive selection
    badSegmentsIndices = false(size(signal));
    prevIdx = [];  % Store the previous index for interpolation
    selecting = true;  % Track whether we are selecting or deselecting

    % Add a button to finalize selection
    uicontrol('Style', 'pushbutton', 'String', 'REMOVE', ...
              'FontWeight', 'bold', 'Position', [20 20 100 40], ...
              'Callback', @finalizeSelection);
    
    % Set the button down function for the figure
    set(f, 'WindowButtonDownFcn', @mouseDown);
    set(f, 'WindowButtonUpFcn', @mouseUp);

    % Mouse down callback function
    function mouseDown(~, ~)
        cp = get(ax, 'CurrentPoint');
        xClick = cp(1,1);
        [~, idx] = min(abs(t - xClick));
        selecting = ~badSegmentsIndices(idx);  % Determine if we are selecting or deselecting
        set(f, 'WindowButtonMotionFcn', @mouseMove);
        mouseMove();  % Capture the initial point as well
    end

    % Mouse move callback function
    function mouseMove(~, ~)
        % Get the current point
        cp = get(ax, 'CurrentPoint');
        xClick = cp(1,1);
        
        % Find the closest time index
        [~, idx] = min(abs(t - xClick));
        
        % If previous index exists, interpolate between the previous and current indices
        if ~isempty(prevIdx)
            if idx > prevIdx
                newIndices = prevIdx:idx;
            else
                newIndices = idx:prevIdx;
            end
        else
            newIndices = idx;
        end
        
        % Update selection state based on the selecting flag
        badSegmentsIndices(newIndices) = selecting;

        % Update previous index
        prevIdx = idx;

        % Update plot to highlight bad segments
        tmpSignal = signal;
        tmpSignal(badSegmentsIndices) = NaN;
        
        % Clear previous plot and update
        cla(ax);
        plot(ax, t, tmpSignal, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
        hold(ax, 'on');
        plot(ax, t(badSegmentsIndices), signal(badSegmentsIndices), '.', ...
            'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1);  % Plot bad segments in orange
        title(ax, 'Select bad segments to remove');
        xlabel(ax, 'Time (minutes)', 'FontWeight', 'bold');
        ylabel(ax, 'Amplitude', 'FontWeight', 'bold');
    end

    % Mouse up callback function
    function mouseUp(~, ~)
        set(f, 'WindowButtonMotionFcn', '');
        prevIdx = [];  % Reset the previous index
    end

    % Finalize selection callback function
    function finalizeSelection(~, ~)
        % Remove bad segments from the signal and time index
        keepIndices = ~badSegmentsIndices;
        cleanSignal = signal(keepIndices);
        % Adjust cleanT to be continuous from 0 to the length of cleanSignal
        cleanT = linspace(t(1), t(end), length(cleanSignal));

        % Close the figure
        delete(f);

        % Display final message
        disp('All bad segments removed and cleaned signal plotted.');
    end

    % Close figure callback function
    function closeFigure(~, ~)
        % Close the figure safely
        delete(f);
        disp('Figure closed without error.');
    end

    % Wait for user interaction before continuing
    uiwait(f);

end
