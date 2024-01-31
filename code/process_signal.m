%% Filter and plot electroencephalography (EGG) signal
clear; close all; clc

mainDir = 'G:\Shared drives\Grants\Granters (Foundations + Funders)\Bial\2022\(000) Yount_Bial_2022\Telly Belly Research';
codeDir = fullfile(mainDir, 'eeg_code');
dataDir = fullfile(mainDir, 'tests');
cd(dataDir)
eeglab; close;

filename = 'test_003.edf';
newFs = 100;  % downsampling to this freq (in Hz)
fc = 0.1;     % cutoff freq for lowpass filter (in Hz)

% CSV file
% tmp = readmatrix(fullfile(dataDir, filename));
% t = tmp(:,1)';
% signal = tmp(:,2)';

% EDF file
EEG = import_edf(fullfile(dataDir,filename));

% Remove bad segments with large artifacts
% pop_eegplot(EEG,1,1,1);
if str2double(filename(8)) == 3
    EEG = eeg_eegrej(EEG, [892 339384;1454636 1542091;1735657 1882719;2160500 2188274;2222102 2244000]);
elseif str2double(filename(8)) == 5
    EEG = eeg_eegrej( EEG, [3370505 4802000]);
% elseif str2double(filename(8)) == 6
%     EEG = eeg_eegrej( EEG, [32751 36500] );
end

% Remove offsets from data rejection

% convert to minutes
t = EEG.times ./ 1000 ./ 60;
signal = EEG.data;

% Interpolate NaNs if any
tf = isnan(signal);
if sum(tf)>0
    warning('%g NaNs detected --> interpolating them',sum(tf))
    ix = 1:length(signal);
    signal(tf) = interp1(ix(~tf),signal(~tf),ix(tf));
end

% Plot raw signal
figure('color','w'); 
subplot(2,1,1)
% envelope(signal,fs*5,'peak')
plot(t,signal,'k-.')
hold on; 

% Sample rate and Nyquist freq (times must be minutes!)
total_time = t(end)-t(1);
nSamples = length(signal);
fs = round(nSamples / (total_time*60)); % convert to seconds as it is in minutes
nf = fs/2;  % Nyquist freq

% Downsample/decimate
fac = fs / newFs; % downsample factor
if fac ~= floor(fac)
    fac = round(fac);
    signal = decimate(signal, fac);
else
    signal = resample(signal, 1, fac);
end
nf = newFs/2;  % update Nyquist
t = t(1:fac:end);
fs = newFs;

% Lowpass filter
Wn = fc/nf;                   % normalized cutoff frequency
[z,p,k] = butter(3,Wn,'low');  % Butterworth filter
[sos,g] = zp2sos(z,p,k);
% freqz(sos, 2^14, fs)
signal = filtfilt(sos,g,signal);
plot(t,signal,'r','linewidth',1)
title(sprintf('Raw EGG: %s',filename(1:end-4))); 
xlabel('Time (min)')
legend('raw','lowpass filtered at 0.1 Hz')
axis tight

% Smooth using a span of 10% of data points
% tic
% sig = smooth(t2,sig,0.1,'rloess');
% toc

% % Highpass filter
% fc = 0.005;                       % cutoff freq
% Wn = fc/nf;                       % normalized cutoff frequency
% [z,p,k] = butter(3,Wn,'high');    % Butterworth filter
% [sos,g] = zp2sos(z,p,k);
% % freqz(sos, 2^14, Fs)
% signal = filtfilt(sos,g,signal);
% hold on
% plot(t, signal,'b','linewidth',1.5)

% Add bar to indicate eating event (test_005)
if str2double(filename(8)) == 5
    y_limits = ylim; 

    % Eating start
    idx = find(round(t)==11, 1);
    line([t(idx) t(idx)], y_limits,'color','black','linewidth',1)
    text(t(idx), y_limits(2), 'Eating start', ...
        'HorizontalAlignment','right','VerticalAlignment','top','Color','black')

    % Eating end
    idx = find(round(t)==16, 1);
    line([t(idx) t(idx)], y_limits,'color','black','linewidth',1)
    text(t(idx), y_limits(2), 'Eating end', ...
        'HorizontalAlignment','right','VerticalAlignment','top','Color','black')

    legend('raw','lowpass filtered at 0.1 Hz', '', '')
end

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

print(gcf, fullfile(dataDir,sprintf('%s.png',filename(1:end-4))),'-dpng','-r300');   % 300 dpi .png

%% Truncate file test_005

figure('color','w')

% Before eating
idx = find(round(t)==11, 1);
t2 = t(1:idx);
signal2 = signal(1:idx);
subplot(3,2,1)
plot(t2,signal2,'r','linewidth',1)
title("Before eating"); xlabel('Time (min)'); axis tight
n = length(signal2); 
times = (0:n-1) / fs;
fmin = 0.005; fmax = 0.1;
freqs = linspace(fmin, fmax, 1000);
[power, f] = plomb(signal2, times, freqs, 'normalized');
f = f * 60; % convert to cpm (cycles per minute)
subplot(3,2,2)
plot(f, power)
title('Lomb-Scargle Periodogram')
xlabel('Frequency (cpm)'); ylabel('Normalized Power'); axis tight

% During eating
idx1 = find(round(t)==11, 1);
idx2 = find(round(t)==16, 1);
t2 = t(idx1:idx2);
signal2 = signal(idx1:idx2);
subplot(3,2,3)
plot(t2,signal2,'r','linewidth',1)
title("After eating"); xlabel('Time (min)'); axis tight
n = length(signal2); 
times = (0:n-1) / fs;
fmin = 0.005; fmax = 0.1;
freqs = linspace(fmin, fmax, 1000);
[power, f] = plomb(signal2, times, freqs, 'normalized');
f = f * 60; % convert to cpm (cycles per minute)
subplot(3,2,4)
plot(f, power)
title('Lomb-Scargle Periodogram')
xlabel('Frequency (cpm)'); ylabel('Normalized Power'); axis tight

% After eating
idx = find(round(t)==16, 1);
t2 = t(idx:end);
signal2 = signal(idx:end);
subplot(3,2,5)
plot(t2,signal2,'r','linewidth',1)
title("After eating"); xlabel('Time (min)'); axis tight
n = length(signal2); 
times = (0:n-1) / fs;
fmin = 0.005; fmax = 0.1;
freqs = linspace(fmin, fmax, 1000);
[power, f] = plomb(signal2, times, freqs, 'normalized');
f = f * 60; % convert to cpm (cycles per minute)
subplot(3,2,6)
plot(f, power)
title('Lomb-Scargle Periodogram')
xlabel('Frequency (cpm)'); ylabel('Normalized Power'); axis tight

print(gcf, fullfile(dataDir,sprintf('%s_truncated.png',filename(1:end-4))),'-dpng','-r300');   % 300 dpi .png


%% EEGLAB approach

% eeglab; close;
% EGG = import_edf(fullfile(dataDir, 'test_004.edf'), 1);
% EGG = pop_resample(EGG,100);
% 
% % Remove large artifacts
% % EGG = eeg_eegrej(EGG, [1 24041;70932 77665;85981 94061;107209 112200]);
% % eegplot(EGG.data,'srate',EGG.srate,'spacing',1, 'winlength',60,'plottitle','Raw EGG');
% 
% % Filter after to deal with large discontinuities
% EGG = pop_eegfiltnew(EGG,'locutoff',0.005);
% EGG = pop_eegfiltnew(EGG,'hicutoff',0.1);
% 
% % Plot raw (without DC drift) signal
% figure('color','w'); 
% subplot(2,1,1);
% plot(EGG.times./1000/60,EGG.data,'linewidth',1.5); xlabel('Time (min)')
% axis tight
% 
% % Lomb-Scargle periodogram
% band = [0.008 0.1];
% nfft = EGG.srate*10;    % 30 s window
% fvec = band(1):1/nfft:band(2);
% [pwr,freqs] = plomb(EGG.data,EGG.times,fvec,"psd"); 
% subplot(2,1,2)
% plot(freqs,pwr); xlabel('Frequency (Hz)'); ylabel('Power'); title('Lomb-Scargle periodogam')
