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

%Peak at 0 Hz due to DC mean component
signal = detrend(signal); %detrend data to remove DC Mean Component (subtract average from each 
                                      %data-point) - Now there is no peak at 0 Hz 

signal_HannWnd = signal.*hanning(L, 'periodic'); %Apply Hann Window
signal_HannWnd = padarray(signal_HannWnd, padding); %zero padding of signal
signaldft_HannWnd = fftshift(fft(signal_HannWnd,nfft)); %FFT of signal
msignaldft = abs(signaldft_HannWnd); %Absolute values of FFT

f = Fs*(-nfft/2:nfft/2-1)/nfft; %Frequency Range

dt = 1/Fs; %Frequency resolution

figure(1)
set(0,'defaultaxesFontName', 'CMU Serif Roman')
set(0,'defaultaxesFontSize', 12)
%set(gcf, 'PaperPositionMode', );
subplot(2,1,1)
plot(((0:(n-1)) * dt),signal)
title('Response');
xlabel('Time (s)'); 
ylabel('Amplitude (degrees)');
subplot(2,1,2)  
plot(f,2*msignaldft);
zoom xon;
zoom(100);
grid on
title('Amplitude Spectrum');
xlabel('Frequency (Hz)'); 
ylabel('Amplitude (dB)');
legend(sprintf( '%s\n%s', 'FFT of Response', 'with Hanning Window'), 'Location', 'NorthEast');
legend('boxoff');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 6];
print -dpng FreqResp.png
end