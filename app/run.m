%% Visualize electrograstrography (EGG) signal - executable program for Garry
%
% Cedric Cannard, June 2024

% clear; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hard-coded parameters %%%%%%%%%%%%%%%%%%%%%%%%
fs          = 100;      % downsample to this freq (in Hz)
lowpass     = 0.1;      % low-pass filter (in Hz)
highpass    = 0.005;    % high-pass filter (in Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % Add path to folder
% tmp = fileparts(which('egg_app'));
% addpath(tmp)

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
    nyquist = 0.5;  % Nyquist frequency
    cutoff = fc * nyquist;
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
tbw = .25*pb;           % transition bandwidth = 25% of pb
cutoff = pb-tbw/2;      % cutoff freq
m = est_filt_order('hamming',fs,tbw);
[b, a] = design_filt(m, cutoff/(fs/2),'high'); % design filter
if m>1000
    fprintf("Filter order (%g) is greater than 1000. Using frequency domain to filter signal faster... \n",m)
    signal = fir_filterdcpadded(b, a, signal', 0, 1)';  % forward/reverse digital filtering in frequency domain

else
    signal = filtfilt(b,a,signal);  % forward/reverse digital filtering in time domain
end

% figure; plot(signal,'k--'); hold on; plot(signal_filt,'b')

% Highpass filter
disp("Performing low-pass filter to remove slow-frequency drifts...")
pb = lowpass;          % pass-band edge
tbw = .25*pb;           % transition bandwidth = 25% of pb
cutoff = pb-tbw/2;      % cutoff freq
m = est_filt_order('hamming',fs,tbw);
[b, a] = design_filt(m, cutoff/(fs/2),'low'); % design filter
if m>1000
    fprintf("Filter order (%g) is greater than 1000. Using frequency domain to filter signal faster... \n",m)
    signal = fir_filterdcpadded(b, a, signal', 0, 1)';  % forward/reverse digital filtering in frequency domain
else
    signal = filtfilt(b,a,signal);  % forward/reverse digital filtering in time domain
end

 
% convert time to minutes
t = t ./ 1000 ./ 60;

% Remove bad segments
[cleanSignal, cleanT] = clean_signal(signal, t);
signal = cleanSignal; t = cleanT;

% Plot raw signal
figure('color','w','NumberTitle', 'off','MenuBar', 'none', 'ToolBar', 'none'); 
subplot(2,1,1)
% envelope(signal,fs*5,'peak')
% plot(t,signal,'k-.')
plot(t,signal,'k','LineWidth',1)
hold on; axis tight; 
title(sprintf("Raw time series - File %s ",fileName(1:end-4)));
xlabel("Time (min)"); ylabel('Amplitude')
% fs = EGG.srate;


% Lomb-Scargle Periodogram
disp("Computing Lomb-scargle periodogram...")
n = length(signal);
times = (0:n-1) / fs;
fmin = 0.005; % min freq in Hz (e.g., 0.005 Hz)
fmax = 0.1;   % max freq in Hz (e.g., 0.1 Hz)
freqs = linspace(fmin, fmax, 1000); % 1000 frequency points
[power, f] = plomb(signal, times, freqs, 'normalized');
f = f * 60; % convert to cpm (cycles per minute)
subplot(2,1,2)
plot(f, power)
title('Lomb-Scargle Periodogram')
xlabel('Frequency (cpm)'); ylabel('Normalized Power')
axis tight
set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

fprintf("Saving plot in the same location: %s \n", fullfile(filePath,sprintf('%s_power-spectrum.png',fileName(1:end-4))))
print(gcf, fullfile(filePath,sprintf('%s_power-spectrum.png',fileName(1:end-4))),'-dpng','-r300');   % 300 dpi .png



%% subfunctions

