%% Visualize electrograstrography (EGG) signal - executable program for Garry
%
% Cedric Cannard, June 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hard-coded parameters %%%%%%%%%%%%%%%%%%%%%%%%
fs          = 100;      % downsample to this freq (in Hz)
lowpass     = 0.1;      % low-pass filter (in Hz)
highpass    = 0.005;    % high-pass filter (in Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select file on computer
[fileName, filePath] = uigetfile({ '*.edf' }, 'Select .edf file');
fullPath = fullfile(filePath, fileName);
if fileName==0
    disp('Aborted.'); return
end

% Load data
[signal, t, sRate, sig_len] = load_data(fullPath);

% Downsample/decimate if sample rate is >100 Hz
if sRate > 100
    fac = sRate / fs;  % Downsampling factor

    % Anti-aliasing filter parameters
    fc = 0.9;  % Cutoff frequency of the anti-aliasing filter (normalized frequency)
    df = 0.2;  % Transition band width of the anti-aliasing filter (normalized frequency)

    % Design the anti-aliasing filter (FIR with Kaiser window)
    cutoff = 0.5*fc; % Nyquist frequency
    n = 100;  % Filter order (adjust as needed)
    b = fir1(n, cutoff, 'low', kaiser(n+1, 5));
    signal = filtfilt(b, 1, signal); % zero-phase filtering

    % Decimate the filtered signal with downsampling factor fac
    if fac ~= floor(fac)
        signal_resamp = decimate(signal, round(fac));
    else
        signal_resamp = resample(signal, 1, fac);
    end

    % Adjust time index accordingly
    t_resamp = linspace(t(1), t(end), length(signal_resamp));

    % Overwrite original variables
    signal = signal_resamp;
    t = t_resamp;
    sRate = fs;  % Update the sample rate to the new rate after resampling
    sig_len = length(signal);  % Update the signal length
end

% Highpass filter
disp("Performing high-pass filter to remove slow-frequency drifts...")
pb = highpass;          % pass-band edge
% tbw = .5*pb;           % transition bandwidth = 50% of pb
tbw = pb;          % transition bandwidth = pb as in pop_eegfiltnew()
cutoff = pb-tbw/2;      % cutoff freq
m = est_filt_order('hamming',fs,tbw);
[b, a] = design_filt(m, cutoff/(fs/2),'high'); % design filter
if m>1000
    fprintf("Filter order (%g) is greater than 1000. Using frequency domain to filter signal faster... \n",m)
    signal = fir_filterdcpadded(b, a, signal', 0, 1)';  % forward/reverse digital filtering in frequency domain

else
    signal = filtfilt(b,a,signal);  % forward/reverse digital filtering in time domain
end

% Lowpass filter
disp("Performing low-pass filter to remove slow-frequency drifts...")
pb = lowpass;          % pass-band edge
% tbw = .5*pb;           % transition bandwidth = 50% of pb
tbw = pb;           % transition bandwidth = pb as in pop_eegfiltnew()
cutoff = pb+tbw/2;      % cutoff freq
m = est_filt_order('hamming',fs,tbw);
[b, a] = design_filt(m, cutoff/(fs/2),'low'); % design filter
if m>1000
    fprintf("Filter order (%g) is greater than 1000. Using frequency domain to filter signal faster... \n",m)
    signal = fir_filterdcpadded(b, a, signal', 0, 1)';  % forward/reverse digital filtering in frequency domain
else
    signal = filtfilt(b,a,signal);  % forward/reverse digital filtering in time domain
end


% Convert time to minutes
t = t ./ 1000 ./ 60;

% Remove bad segments
[cleanSignal, cleanT] = clean_signal(signal, t);
signal = cleanSignal; t = cleanT;

% Plot raw signal
figure('color','w','NumberTitle', 'off','MenuBar', 'none', 'ToolBar', 'none');
subplot(2,1,1)
plot(t,signal,'k','LineWidth',1)
hold on; axis tight;
title(sprintf("Raw time series - File %s ",fileName(1:end-4)));
xlabel("Time (min)"); ylabel('Amplitude')

% Lomb-Scargle Periodogram
disp("Computing Lomb-scargle periodogram...")
fmin = highpass;    % min freq in Hz
fmax = lowpass;     % max freq in Hz
freqs = linspace(fmin, fmax, 1000); % 1000 frequency points
[power, f] = plomb(signal, t*60, freqs, 'normalized');
f = f * 60; % convert to cpm (cycles per minute)
subplot(2,1,2)
plot(f, power)
title('Lomb-Scargle Periodogram')
xlabel('Frequency (cpm)'); ylabel('Normalized Power')
axis tight
set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

fprintf("Saving plot in the same location: %s \n", fullfile(filePath,sprintf('%s_power-spectrum.png',fileName(1:end-4))))
print(gcf, fullfile(filePath,sprintf('%s_power-spectrum.png',fileName(1:end-4))),'-dpng','-r300');   % 300 dpi .png


%% SUBFUNCTIONS

% LOAD_DATA Load EEG data from a selected .edf file and extract relevant information.
%
% Usage:
%   [signal, t, fs, len] = load_data()
%
% Inputs:
%   None. The function uses a file selection dialog to choose a .edf file.
%
% Outputs:
%   signal - A matrix representing the EEG signal data (channels x samples).
%   t      - A vector representing the time points corresponding to the signal.
%   fs     - The sampling rate of the EEG data.
%   len    - The length of the EEG signal (number of samples).
%
% Example:
%   % Load EEG data from a selected .edf file
%   [signal, t, fs, len] = load_data();
%
%   % Plot the first 10 seconds of the signal
%   figure;
%   plot(t(1:fs*10), signal(1,1:fs*10));
%   xlabel('Time (seconds)');
%   ylabel('Amplitude');
%   title('EEG Signal (First 10 Seconds)');
%
% Copyright (C), Cedric Cannard, June 2024

function [signal, t, fs, len] = load_data(fullPath)

% Load file and extract info
fprintf("Loading selected file %s... \n", fullPath);
[tmp_signal, annot] = edfread(fullPath, 'TimeOutputType', 'datetime');
info = edfinfo(fullPath);
annot = timetable2table(annot,'ConvertRowTimes',true);

% Timestamps in datetime format
edfTime = timetable2table(tmp_signal,'ConvertRowTimes',true);
edfTime = datetime(table2array(edfTime(:,1)), 'Format', 'HH:mm:ss:SSS');

% Sampling rate
sPerCell = mode(seconds(diff(edfTime)));
if sPerCell == 1
    fs = info.NumSamples(1);
else
    fs = info.NumSamples(1)/sPerCell;
end

% Signal
tmp_signal = table2array(tmp_signal)';

% For test EEG dataset: keep only first channel
if size(tmp_signal,1)>1
    tmp_signal = tmp_signal(1,:);
end

% Convert from cells to one double-precision array
signal = [];
for iChan = 1:size(tmp_signal,1)
    sample = 1;
    for iCell = 1:size(tmp_signal,2)
        cellData = tmp_signal{iChan,iCell};
        if sPerCell == 1     %data with correct sample rate at import
            signal(iChan, sample:sample+fs-1) = cellData;
            sample = sample + fs;
        else
            % data with incorrect sample rate at import (e.g. RKS05 or RKS09)
            for iSec = 1:sPerCell
                if iSec == 1
                    signal(iChan, sample:sample+fs-1) = cellData(iSec:iSec*fs);
                    sample = sample + fs;
                else
                    signal(iChan, sample:sample+fs-1) = cellData(((iSec-1)*fs)+1 : (iSec)*fs);
                    sample = sample + fs;
                end
            end
        end
    end
end
clear tmp_signal


% Time index
len = size(signal,2);
% len = sample;
t = ( (0:len-1) / fs ) .* 1000;  % in ms
% edfTime = edfTime - edfTime(1);
% t = round(linspace(milliseconds(edfTime(1)), milliseconds(edfTime(end)), len),5);
% times(1:5)
% t(1:5)


% print the total duration (in minutes + seconds)
len_sec = t(end)/1000;
fprintf('Total duration: %g min  \n', len_sec);
fprintf('Total duration: %g min %g sec. \n', floor(len_sec/60), round(mod(len_sec,60)));

end



% Estimate windowed sinc FIR filter order depending on window type and
% requested transition band width
% 
% Usage:
%   [m, dev] = firwsord(wtype, fs, df);
%   m = firwsord('kaiser', fs, df, dev);
%
% Inputs:
%   wtype - char array window type. 'rectangular', 'hann', 'hamming',
%           'blackman', or 'kaiser'
%   fs    - scalar sampling frequency
%   df    - scalar requested transition band width
%   dev   - scalar maximum passband deviation/ripple (Kaiser window
%           only)
%
% Output:
%   m     - scalar estimated filter order
%   dev   - scalar maximum passband deviation/ripple
%
% References:
%   [1] Smith, S. W. (1999). The scientist and engineer's guide to
%       digital signal processing (2nd ed.). San Diego, CA: California
%       Technical Publishing.
%   [2] Proakis, J. G., & Manolakis, D. G. (1996). Digital Signal
%       Processing: Principles, Algorithms, and Applications (3rd ed.).
%       Englewood Cliffs, NJ: Prentice-Hall
%   [3] Ifeachor E. C., & Jervis B. W. (1993). Digital Signal
%       Processing: A Practical Approach. Wokingham, UK: Addison-Wesley

function m = est_filt_order(wintype,fs,df,dev)

winTypeArray = {'rectangular', 'hann', 'hamming', 'blackman', 'kaiser'};
winDfArray = [0.9 3.1 3.3 5.5];
% winDevArray = [0.089 0.0063 0.0022 0.0002];

wintype = find(strcmp(wintype, winTypeArray));
df = df / fs; % Normalize transition bandwidth

if strcmpi(wintype,'kaiser') % Kaiser window
    devdb = -20*log10(dev);
    m = 1 + (devdb - 8) / (2.285 * 2 * pi * df);
else
    m = winDfArray(wintype) / df;
    % dev = winDevArray(wintype);
end

% Make filter order even (FIR type I)
m = ceil(m / 2) * 2;

end

% Andreas Widmann function
function data = fir_filterdcpadded(b, a, data, causal, usefftfilt)

% Defaults
if nargin < 5 || isempty(usefftfilt)
    usefftfilt = 0;
end
if nargin < 4 || isempty(causal)
    causal = 0;
end

% Check arguments
if nargin < 3
    error('Not enough input arguments.');
end

% Is FIR?
if ~isscalar(a) || a ~= 1
    error('Not a FIR filter. onepass-zerophase and onepass-minphase filtering is available for FIR filters only.')
end

% Group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;

% Filter symmetry
isSym = all(b(1:groupDelay) == b(end:-1:groupDelay + 2));
isAntisym = all([b(1:groupDelay) == -b(end:-1:groupDelay + 2) b(groupDelay + 1) == 0]);
if causal == 0 && ~(isSym || isAntisym)
    error('Filter is not anti-/symmetric. For onepass-zerophase filtering the filter must be anti-/symmetric.')
end

% Padding
if causal
    startPad = repmat(data(1, :), [2 * groupDelay 1]);
    endPad = [];
else
    startPad = repmat(data(1, :), [groupDelay 1]);
    endPad = repmat(data(end, :), [groupDelay 1]);
end

if usefftfilt
    data = fftfilt(double(b), double([startPad; data; endPad]));
else
    data = filter(double(b), 1, double([startPad; data; endPad])); % Pad and filter with double precision
end

% Remove padded data
data = data(2 * groupDelay + 1:end, :);

end



% Design windowed sinc type I linear phase FIR filter
%
% Usage:
%   b = design_filt(m, f, t, w);
%
% Inputs:
%   m - filter order (mandatory even)
%   f - vector or scalar of cutoff frequency/ies (-6 dB;
%       pi rad / sample)
%
% Optional inputs:
%   w - vector of length m + 1 defining window {default hamming}
%   t - 'high' for highpass, 'stop' for bandstop filter {default low-/
%       bandpass}
%
% Output:
%   b - filter coefficients
%
% References:
%   Smith, S. W. (1999). The scientist and engineer's guide to digital
%   signal processing (2nd ed.). San Diego, CA: California Technical
%   Publishing.

function [b, a] = design_filt(m, f, t, w)

a = 1;

if nargin < 2
    error('Not enough input arguments');
end
if length(m) > 1 || ~isnumeric(m) || ~isreal(m) || mod(m, 2) ~= 0 || m < 2
    error('Filter order must be a real, even, positive integer.');
end
f = f / 2;
if any(f <= 0) || any(f >= 0.5)
    error('Frequencies must fall in range between 0 and 1.');
end
if nargin < 3 || isempty(t)
    t = '';
end
if nargin < 4 || isempty(w)
    if ~isempty(t) && ~ischar(t)
        w = t;
        t = '';
    else
        w = windows('hamming', (m + 1));
    end
end
w = w(:)'; % Make window row vector

b = fkernel(m, f(1), w);

if isscalar(f) && strcmpi(t, 'high')
    b = fspecinv(b);
end

if length(f) == 2
    b = b + fspecinv(fkernel(m, f(2), w));
    if isempty(t) || ~strcmpi(t, 'stop')
        b = fspecinv(b);
    end
end
end


% Compute filter kernel
function b = fkernel(m, f, w)
m = -m / 2 : m / 2;
b(m == 0) = 2 * pi * f; % No division by zero
b(m ~= 0) = sin(2 * pi * f * m(m ~= 0)) ./ m(m ~= 0); % Sinc
b = b .* w; % Window
b = b / sum(b); % Normalization to unity gain at DC
end


% Spectral inversion
function b = fspecinv(b)
b = -b;
b(1, (length(b) - 1) / 2 + 1) = b(1, (length(b) - 1) / 2 + 1) + 1;
end

% windows() - Symmetric window functions
%
% Usage:
%   >> h = windows(t, m);
%   >> h = windows(t, m, a);
%
% Inputs:
%   t - char array 'rectangular', 'bartlett', 'hann', 'hamming',
%       'blackman', 'blackmanharris', 'kaiser', or 'tukey'
%   m - scalar window length
%
% Optional inputs:
%   a - scalar or vector with window parameter(s)
%
% Output:
%   w - column vector window
%
% Author: Andreas Widmann, University of Leipzig, 2014

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2014 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% $Id$

function w = windows(t, m, a)

if nargin < 2 || isempty(t) || isempty(m)
    error('Not enough input arguments.');
end

% Check window length
if m ~= round(m)
    m = round(m);
    warning('firfilt:nonIntegerWindowLength', 'Non-integer window length. Rounding to integer.')
end
if m < 1
    error('Invalid window length.')
end

% Length 1
if m == 1
    w = 1;
    return;
end

% Even/odd?
isOddLength = mod(m, 2);
if isOddLength
    x = (0:(m - 1) / 2)' / (m - 1);
else
    x = (0:m / 2 - 1)' / (m - 1);
end

switch t
    case 'rectangular'
        w = ones(length(x), 1);
    case 'bartlett'
        w = 2 * x;
    case 'hann'
        a = 0.5;
        w = a - (1 - a) * cos(2 * pi * x);
    case 'hamming'
        a = 0.54;
        w = a - (1 - a) * cos(2 * pi * x);
    case 'blackman'
        a = [0.42 0.5 0.08 0];
        w = a(1) - a(2) * cos (2 * pi * x) + a(3) * cos(4 * pi * x) - a(4) * cos(6 * pi * x);
    case 'blackmanharris'
        a = [0.35875 0.48829 0.14128 0.01168];
        w = a(1) - a(2) * cos (2 * pi * x) + a(3) * cos(4 * pi * x) - a(4) * cos(6 * pi * x);
    case 'kaiser'
        if nargin < 3 || isempty(a)
            a = 0.5;
        end
        w = besseli(0, a * sqrt(1 - (2 * x - 1).^2)) / besseli(0, a);
    case 'tukey'
        if nargin < 3 || isempty(a)
            a = 0.5;
        end
        if a <= 0 % Rectangular
            w = ones(length(x), 1);
        elseif a >= 1 % Hann
            w = 0.5 - (1 - 0.5) * cos(2 * pi * x);
        else
            mTaper = floor((m - 1) * a / 2) + 1;
            xTaper = 2 * (0:mTaper - 1)' / (a * (m - 1)) - 1;
            w = [0.5 * (1 + cos(pi * xTaper)); ones(length(x) - mTaper, 1)];
        end
    otherwise
        error('Unkown window type')
end

% Make symmetric
if isOddLength
    w = [w; w(end - 1:-1:1)];
else
    w = [w; w(end:-1:1)];
end

end



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
