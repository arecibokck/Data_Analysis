function [freqs,Y] = fastFourierT(timeVector,y)
% This function returns both properly defined frequency vector and fourier
% transform of a input time and signal vector for real valued signals. 
% 
% call:
% [freqs,Y] = fastFourierT(timeVector,y)
%
% -----Output-----
% freqs         frequency vector
% Y             fft of y
%
% -----Input------
% timeVector    equally spaced vector with times
% y             real space signal
%
% see also custom function: inverseFFT

% normal fft

Y = fft(y)/numel(y);
Y = Y(1:round(end/2)+1);
Y(2:end-1) = 2*Y(2:end-1);
freqs = 1/mean(diff(timeVector)) * (0:round(numel(y)/2)) /round(numel(y));


end