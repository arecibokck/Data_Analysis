function measuredPSD = makeDSH_PSD(frequencies,parameters,varargin)
% This function simulates the power spectral density (PSD) measured in a
% delayed self-heterodyning (DSH) experiment.
% 
% call:
% measuredPSD = makeDSH_PSD(frequencies,parameters)
%
% ... to modify internal variables an additional variable input in the form
% of makeDSH_PSD(...,'variableName',value) can be made.
%
% ---------- Output ---------
% measuredPSD           simulated values of the power spectral density for
%                       the given input frequency sampling and parameters
% 
% ---------- Input ----------
% frequencies           equidistant real valued vector of frequency
%                       sampling points
% parameters            vector including the parameters of the simulation
%                       parameters(1):  beat frequency of the heterodyne
%                                       experiment
%                       parameters(2):  length of the fiber loop (used to
%                                       calculate the time delay)
%                       parameters(3):  amplitude of the beat signal
%                       parameters(4):  amplitude of the Gaussian (white)
%                                       noise contribution
%                       parameters(5):  amplitude of the 1/f (flicker)
%                                       noise
%                       parameters(6):  amplitued of the 1/f^2 (random
%                                       walk) noise
%                       note: the frequency dependence of the noise refers
%                       to the dependence of the frequency noise spectral
%                       density P_\nu\nu 
%
% varargin: some internal variables can be modified via additional input
% pairs of the form 'variable name',value - in the following, the default
% values are shown: 
% 'c',299792458         % the speed of light ... not sure why you would
%                       % like to change that 
% 'ts', ...             % time sampling constant (time difference between
%                       % simulated signal trace points), calculated as
%                       % 1/(5*parameters(1)) to ensure fine enough
%                       % sampling of the beat signal 
% 'numSamp', 1e7        % number of points in time for the simulated signal
%                       % traces - trade-off between simulation time and
%                       % smoothness of the simulated spectrum - stronger
%                       % noise might also indicate the necessity to
%                       % increase this
% 'delayt', ...         % delay time between the two arms of the heterodyne
%                       % experiment (calculated as parameters(2)/c)
% 'RBW',2e3             % resolution bandwidth of the electrical spectrum
%                       % analyzer used in the experiment (Gaussian
%                       % response function assumed) 
% convFactor,1.5116     % conversion factor between the RBW and the width
%                       % of the Gaussian filter performed through the
%                       % convolution with the response function

%% just for testing ... 
% if set to true - execute as script and see comparison plots at end 
testing = false;
if testing
    clearvars -except testing
%     close all

    parameters = [80e6,6.5e3,1,1e3,0,0];
    frequencies = linspace(parameters(1)-300e3, ...
        parameters(1)+300e3,1e3+1);
end

%% read-in simulation parameters
deltaf = parameters(1);                             % beat frequency
looplength = parameters(2);                         % fiber loop length
signalAmplitude = parameters(3);                    % amplitude of the beat signal
noiseAmplitudes = parameters(4:6);                  % noise amplitudes

%% other parameters
c = 299792458;              % speed of light
ts = 1/(5*deltaf);          % time sampling constant
numSamp = 1e7;              % number of time samples
delayt = looplength/c;      % delay time due to loop
RBW = 2e3;                  % resolution bandwidth of the spectrum analyzer
convFactor = 1.4;           % conversion factor between RBW and Gaussian width

%% overwrite defaults - if specified
if ~testing
    for Index = 1:numel(varargin)/2
        eval([varargin{2*Index-1} '= varargin{2*Index};']);
    end
end

%% create incident psd signal
% phase noise:
[timeVector,phaseNoise] = frequencyNoise2(ts,numSamp,noiseAmplitudes);

% loop delay:
[~,delayIndexShift] = min(abs(timeVector-delayt));

% time vector of the signals:
timeVectorSig = timeVector(1:end-delayIndexShift);

% simulated signal on the photodiode:
signal = signalAmplitude*cos(2*pi*deltaf*timeVectorSig ...
    + phaseNoise(delayIndexShift+1:end) - phaseNoise(1:end-delayIndexShift));

% welch estimation of the incident power spectral density - sampling cutoff
% at frequency spacing
cffreq = mean(diff(frequencies));
[freqs,signalSpecDens] = ownWelch(timeVectorSig,signal,'cffreq',cffreq);

% cut to relevant spectrum part ... the frequency window of the queried spectral points
freqWindow = [frequencies(1)-cffreq frequencies(end)+cffreq];
if any((freqs>freqWindow(1))&(freqs<freqWindow(2)))
    signalSpecDens = signalSpecDens((freqs>freqWindow(1))&(freqs<freqWindow(2)));
    freqs = freqs((freqs>freqWindow(1))&(freqs<freqWindow(2)));
end

%% include the resolution bandwidth of the spectrum analyzer
% frequency stepping:
step = mean(diff(freqs));

% window of response function:
DeltaRSFs = -2*RBW:step:2*RBW;

% response function definition
RSF = @(nu) 1./(sqrt(pi)*RBW/convFactor).*exp(-nu.^2./(RBW/convFactor).^2);

% convoluted PSD and response function
signalSpecDensConv = conv(signalSpecDens,RSF(DeltaRSFs),'same').*step;

%% break the simulated spectrum down to the frequency sampling of the input
measuredPSD = interp1(freqs,signalSpecDensConv,frequencies);

%% just for testing ...
if testing
%     figure
    hold on
%     plot(freqs,10*log10(signalSpecDens)); 
%     plot(freqs,10*log10(signalSpecDensConv));
    plot(frequencies,10*log10(measuredPSD)); 
    hold off
    set(gcf, 'WindowStyle', 'Docked');
end

end
