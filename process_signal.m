%% Filter and plot EGG signal
clear; %close all; clc

tmp = readmatrix('C:\Users\Cedric\Downloads\test_EGG_02.csv');
t = tmp(:,1)';
signal = tmp(:,2)';

% Interpolate NaNs if any
tf = isnan(signal);
if sum(tf)>0
    warning('%g NaNs detected --> interpolating them',sum(tf))
    ix = 1:length(signal);
    signal(tf) = interp1(ix(~tf),signal(~tf),ix(tf));
end

% Remove DC drift
ft = fft(signal);
ft(1) = 0;              % zero out the DC component
signal = ifft(ft);      % Inverse transform back to time domain.
if abs(mean(real(signal))) > .0005
    warning('Data mean should be close to 0, DC drift removal must have failed.')
end

% Normalize signal
% signal = zscore(signal);

% sample rate and Nyquist freq
total_time = t(end)-t(1);
nSamples = length(signal);
Fs = round(nSamples / (total_time*60)); %convert to seconds as it is in minutes
nf = Fs/2;  % Nyquist freq

% figure('color','w'); envelope(signal,Fs*5,'peak')

% Downsample/decimate
newFs = 62.5;  % in Hz
fac = Fs / newFs; % downsample factor
if fac ~= floor(fac)
    fac = round(fac);
    sig = decimate(signal, fac);
else
    sig = resample(signal, 1, fac);
end
nf = newFs/2;  % update Nyquist
t2 = t(1:fac:end);

% Lowpass filter
fc = 0.1;                     % cutoff freq
Wn = fc/nf;                   % normalized cutoff frequency
[z,p,k] = butter(12,Wn,'low');  % Butterworth filter
[sos,g] = zp2sos(z,p,k);
% freqz(sos, 2^14, Fs)
sig = filtfilt(sos,g,sig);

% Smooth using a span of 10% of data points
% tic
% sig = smooth(t2,sig,0.1,'rloess');
% toc

% Plot
figure('color','w'); 
subplot(2,1,1)
plot(t,signal,'k-.')
hold on; plot(t2,sig,'r','linewidth',1.5)

% Highpass filter
fc = 0.01;                     % cutoff freq
Wn = fc/nf;                   % normalized cutoff frequency
[z,p,k] = butter(12,Wn,'high');  % Butterworth filter
[sos,g] = zp2sos(z,p,k);
% freqz(sos, 2^14, Fs)
sig = filtfilt(sos,g,sig);

% Add high-passed signal to plot
hold on
plot(t2, sig,'b','linewidth',1.5)
title('EGG signal'); xlabel('Time (min)')
legend('raw','lowpass 0.1 hz', 'highpass 0.01 hz')
axis tight

% Remove 3 first minutes of testing
idx = t2<4;
t2(idx) =[]; 
sig(idx) = [];

subplot(2,1,2)
plot(t2, sig,'b','linewidth',1.5)
title('Cleaner EGG signal (filtered)'); xlabel('Time (min)')
axis tight

% % Lomb-Scargle periodogram
% band = [0.01 0.1];
% nfft = 1024;    % use this instead? 2^nextpow2(length(t2))
% fvec = band(1):1/nfft:band(2);
% [pwr,freqs] = plomb(sig,t2,fvec,"psd"); 
% subplot(3,1,3)
% plot(freqs,pwr)
% xlabel('Frequency (Hz)'); ylabel('Power'); title('Lomb-Scargle periodogam')

gong
