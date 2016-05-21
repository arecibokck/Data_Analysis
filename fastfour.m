% This function outputs frequency vs. amplitude plot of an input given time series
% The function requires the following inputs
% signal=times series data
% sam=frequency of sampling

function [fouranal]=fastfour(response,sam);
signal = csvread(response);
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

signal_HannWnd = signal.*hanning(L, 'periodic'); %Apply Hann Window
signal_HannWnd = padarray(signal_HannWnd, padding); %zero padding of signal
signaldft_HannWnd = fftshift(fft(signal_HannWnd,nfft)); %FFT of signal
msignaldft = abs(signaldft_HannWnd); %Absolute values of FFT

f = Fs*(-nfft/2:nfft/2-1)/nfft; %Frequency Range

dt = 1/Fs; %Frequency resolution

figure(1),
subplot(2,1,1)
plot(((0:(n-1)) * dt),signal)
title('Response');
xlabel('Time (s)'); 
ylabel('Signal(t)');
subplot(2,1,2)  
plot(f,2*msignaldft);
zoom xon;
zoom(150);
title('Amplitude Spectrum with Hann Window');
xlabel('Frequency (Hz)'); 
ylabel('|Signal(f)|');
end