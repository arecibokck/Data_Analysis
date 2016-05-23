% This function outputs frequency vs. amplitude plot of an input given time series
% The function requires the following inputs
% signal=times series data
% sam=frequency of sampling

function [fouranal]=fastfour(stimulus,response,sam)

%Peak at 0 Hz due to DC mean component
signal_s = detrend(csvread(stimulus)); %detrend data to remove DC Mean Component (subtract average from each 
signal_r = detrend(csvread(response)); %data-point) - Now there is no peak at 0 Hz 

Fs = sam;

if (isunix) %# Linux, Mac
    [status, result] = system( ['wc -l ', response] );
    n = str2num(result);

elseif (ispc) %# Windows
    n = str2num( perl('countlines.pl', response) );

else
    error('...');

end

L = n; % Window Length of FFT    
nfft = 2^nextpow2(L); % Transform length
padding = (nfft - L)/2; %zero padding length on either side of signal

msignalsdft = Signal_Mod(signal_s, L, nfft, padding); 
msignalrdft = Signal_Mod(signal_r, L, nfft, padding); 

f = Fs*(-nfft/2:nfft/2-1)/nfft; %Frequency Range

dt = 1/Fs; %Frequency resolution

figure(1)
set(0,'defaultaxesFontName', 'CMU Serif Roman')
set(0,'defaultaxesFontSize', 12)
subplot(2,2,1)
plot(((0:(n-1)) * dt),signal_s)
title('Stimulus');
xlabel('Time (s)'); 
ylabel('Amplitude (degrees)');
subplot(2,2,3)  
plot(f,2*msignalsdft);
zoom xon;
zoom(30);
grid on
title('FFT of Stimulus');
xlabel('Frequency (Hz)'); 
ylabel('Amplitude (dB)');
subplot(2,2,2)
plot(((0:(n-1)) * dt),signal_r)
title('Response');
xlabel('Time (s)'); 
ylabel('Amplitude (degrees)');
subplot(2,2,4)  
plot(f,2*msignalrdft);
zoom xon;
zoom(30);
grid on
title('FFT of Response');
xlabel('Frequency (Hz)'); 
ylabel('Amplitude (dB)');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 6];
print -dpng FourSpec.png

[y_a, lag_a] = xcorr(signal_s, 'coeff')
%[~,A] = max(abs(y_a));
%lagDiff_a = lag_a(A)
%timeDiff_a = lagDiff_a/Fs
t_a = (-length(y_a)/2:length(y_a)/2-1)/Fs

[y_c, lag_c] = xcorr(signal_s, signal_r, 'coeff')
%[~,C] = max(abs(y_c));
%lagDiff_c = lag_c(C)
%timeDiff_c = lagDiff_c/Fs
t_c = (-length(y_c)/2:length(y_c)/2-1)/Fs

figure(2)
set(0,'defaultaxesFontName', 'CMU Serif Roman')
set(0,'defaultaxesFontSize', 12)
plot(t_a, y_a, 'b.');
hold on;
plot(t_c, y_c, 'r.');
title('Correlation');
xlabel('Lag (s)'); 
ylabel('CF');
legend('AutoCorr', 'CrossCorr', 'Location', 'NorthEast');

end

function [signaldft_HannWnd] = Signal_Mod(signal, L, nfft, padding)

signal_HannWnd = signal.*hanning(L, 'periodic'); %Apply Hann Window
signal_HannWnd = padarray(signal_HannWnd, padding); %zero padding of signal
signaldft_HannWnd = abs(fftshift(fft(signal_HannWnd,nfft))); %Absolute values of FFT of signal
return

end