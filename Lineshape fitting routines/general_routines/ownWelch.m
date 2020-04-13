function [freqs,powerSpectrum] = ownWelch(timeVector,signal,varargin)
% This function performs a Welch estimation of the power spectral density
% of a given input signal (using averaging on the estimated power spectral
% densities of overlapping signal segments). 
% 
% call:
% [freqs,powerSpectrum] = ownWelch(timeVector,signal,varargin)
%
% ... to modify internal variables an additional variable input in the form
% of ownWelch(...,'variableName',value) can be made.
%
% ---------- Output ----------
% freqs             vector of frequency values for the spectrum
% powerSpectrum     estimated power spectral density P_ss for an input
%                   signal s (magnitude squared of the fourier transform of
%                   the signal or fourier transform autocorrelated signal
%                   respectively) 
% 
% ---------- Input -----------
% timeVector        equally spaced vector with the corresponding times of
%                   the signal values 
% signal            real vector with the values of the incident signal s
%
% varargin: some internal variables can be modified via additional input
% pairs of the form 'variable name',value - in the following, the default
% values are shown: 
% 'intervalShift',0.5           % overlap of sample intervals used to
%                               % evaluate the DFT/FFT over which the PSD
%                               % is averaged (50% overlap gives equal
%                               % weighting for data point in the Hann
%                               % windowed samples)  
% 'windowFunction','Hann'       % windowing shape of the samples of the
%                               % signal, apart from 'Hann' also 'Gaussian'
%                               % and 'no' are available (Gaussian window:
%                               % falls to 1/e^4 at edges, no: no window =
%                               % rectangular window of full width)
% 'padding','on'                % to each interval sample a zero padding is
%                               % applied before DFT/FFT to interpolate the
%                               % frequency space ... can be set to 'off'
%                               % (= no padding)
% 'padSize',2                   % size of the zero pad compared to the
%                               % sample interval length (2 means two times
%                               % the sample length of zeros attached) 
% 'cffreq', ...                 % positive real-valued number giving the
%                               % cut-off frequency of the mathematically
%                               % reliable frequency space sampling
%                               % (basically the separation of the points
%                               % in frequency space - finer point sampling
%                               % is due to interpolation and does not
%                               % contain any more information)
%                               % ... the default value of cffreq is
%                               % computed as 10/(time length of the input
%                               % time vector), but can be overwritten

%% parameter settings
cffreq = 10/(timeVector(end)-timeVector(1));      
                        % sample frequency 
                        % - no noise information finer than that

intervalShift = 0.5;    % overlap of intervals - preferred relative overlap 
                        % of intevals (might be adjusted towards smaller 
                        % overlap depending on time series length)
                        
windowFunction = 'Hann';    % 'Hann' window, 'Gaussian' window and 'no' 
                            % window available

padding = 'on';         % zero padding of intervals
padSize = 2;            % size of the zero padding in units of the interval size

%% overwrite defaults - if specified
for Index = 1:numel(varargin)/2
    eval([varargin{2*Index-1} '= varargin{2*Index};']);
end

%% find fft intervals
numberOftimeSteps = numel(timeVector);
timeSize = timeVector(numberOftimeSteps)-timeVector(1);
timeSteps = timeSize/(numberOftimeSteps-1);

minTimeLength = 1./cffreq;
intervalLength = floor(minTimeLength/timeSteps);
numberOfIndependentIntervals = timeSize/minTimeLength;

if numberOfIndependentIntervals<(1 + (1-intervalShift))
    warning('Incident in Welch estimation: Current lower cut-off frequency not met.')
    cffreq = 1/(timeVector(end)-timeVector(1)); 
    disp(['Current cut-off: ' num2str(cffreq) '. No averaging.']);
    intervals = [1 numberOftimeSteps];
    numberOfIntervals = 1;
    intervalLength = numberOftimeSteps;
else
    numberOfIntervals = floor((numberOftimeSteps-intervalShift*intervalLength)/(intervalLength*(1-intervalShift)));
    intervalLength = floor(numberOftimeSteps/(numberOfIntervals*(1-intervalShift)+intervalShift));
    overlap = ceil(intervalShift*intervalLength);
    
    intervals = [(((1:numberOfIntervals)'-1)*(intervalLength-overlap)+1), ...
        ((1:numberOfIntervals)'*(intervalLength-overlap)+overlap)];
end

%% do separate ffts
switch windowFunction
    case 'Gaussian'
        sigma = intervalLength/4;
        window = intervalLength*exp(-(((1:intervalLength)-intervalLength/2)/sigma).^2)/(sqrt(pi)*sigma);
    case 'Hann'
        window = 2*cos(pi*((1:intervalLength)-intervalLength/2)/intervalLength).^2;
    case 'no'
        window = ones(1,intervalLength);
    otherwise
        warning('Windowing not properly specified, taking no window.')
        window = ones(1,intervalLength);
end

accumulatedFFTpow = 0;

for Index = 1:numberOfIntervals
    timeSelection = timeVector(intervals(Index,1):intervals(Index,2));
    signalSelection = signal(intervals(Index,1):intervals(Index,2));
    signalSelection = signalSelection.*window;
    
    if strcmp(padding,'on')
        zeroPadSize = round(padSize*intervalLength);
        timeSelection = [timeSelection,(timeSelection(end)+(1:zeroPadSize)*timeSteps)];
        signalSelection = [signalSelection,zeros(1,zeroPadSize)];
    end
        
    [freqs,currentFFT] = fastFourierT(timeSelection,signalSelection);
    accumulatedFFTpow = accumulatedFFTpow + abs(currentFFT).^2;
end

powerSpectrum = accumulatedFFTpow/numberOfIntervals;

end