clear all 
close all

% This scripts reads-in a DSH data set and tries a simulation based fit of
% the signal. To ensure a good convergence, the initial guess parameters
% need to be sufficiently close to the final values. It is therefore
% srongly recommended to adjust the intial guess until its visual
% appearance matches the data. Read thoroughly through the comments - I
% know there are a lot, but they are intended to be helpful ...
% 
% The concept and routines used here are widely adapted from:
% W. Ma et al, "Laser frequency noise characterization by self-heterodyne
% with both long and short delay," Appl. Opt. 58, 3555-3563 (2019)
% https://doi.org/10.1364/AO.58.003555
% 
% Following sub-functions are required:
% makeDSH_PSD.m
% frequencyNoise2.m
% ownWelch.m
% fastFourierT.m
% filterGeneratedNoise2.m 
% extractLinewidths.m
% 
% To view the performance on the example files, the respective csv files
% need to be present (14.1.csv or 23.1.csv respectively) and the
% corresponding parts in the script have to be uncommented/ commented (see
% sections: "read data" and "define fit function and find initial guess") 

%% delay line and expereimental parameters
% Following specifies the delay of the signal. The looplength is the
% approximate value of the fiber length and the added path the
% correspondingly added optical path (including the refractive index of the
% fiber). These parameters have to be adjusted until the periodicity of the
% simulated pattern matches the experimentally observed. 
looplength = 4.3e3;
addedpath = 1.47*looplength;

c = 299792458;      % speed of light for later use

% Next the used beat frequency has to be specified (so the modulation
% frequency through the AOM in one arm).  
beatFreq = 80e6;

%% other parameters
% These are a number of parameters that are important for the simulation of
% the noisy sample and mostly deal with the simulated duration and sampling
% and a normalization of the noise that has to be adjusted accordingly.
% 
% The number of samples is the size of the time vectors that are used to
% generate signals or noise. If the noise amplitudes are not normalized to
% it, changes of the initial guess will be needed as well! So keep it
% consistent throughout the whole used code and subroutines.
% Also note that much smaller sample numbers are not useful, since this
% increases the amount of noise on the simulated data making it harder to
% fit it (each evaluation of the simulation even under the same conditions
% will then give slightly other spectra). 
numSamp = 1e7;

% The signal sampling tells the used amount of sampling points during an
% oscillation of the beat signal. Basically, it gives the spacing of the
% time vector normalized to the expected signal period. Larger numbers mean
% finer sampling, which comes at the cost of smaller full times and noisier
% spectra, if 'numSamp' is not adjusted accordingly (which would increase
% the simulation time again). Too small numbers (<4) will get you in
% trouble with actually sampling the desired signal. The chosen pair seemed
% to be kind of a sweet spot for the example files .
signalSampling = 5; 

% Thanks to this factor the spectra should not change for a given set of
% noise amplitudes, if the signal sampling or number of samples are
% changed. If you do not want this, set it to 1. 
noiseAmplitudeNormalizationFactor = sqrt(signalSampling*beatFreq/numSamp);

%% read data
% In the first step the data is read into MATLAB and plotted for
% comparison. 

% examples - uncomment one or the other
% ----------------------------------------------------
% example file demoSpectrum_14.1.csv
fileName = 'demoSpectrum_14.1.csv';
% ----------------------------------------------------
% example file demoSpectrum_23.1.csv
% fileName = 'demoSpectrum_23.1.csv';
% ----------------------------------------------------

% the entries of the data start in row 44, first column being the frequency
% second column the power in dBm 
data = csvread(fileName,44,0);

% rename the parts for convenience
freqs = data(:,1);
powerdbm = data(:,2);
power = 10.^(powerdbm/10);

% plot the measured data - I docked the figure for convenience on my
% limited screen size - feel free to adjust anything.
fh = figure; hold on; 
set(fh, 'WindowStyle', 'Docked');
plot(freqs,powerdbm,'.'); 
xlabel('Beat frequency [Hz]');
ylabel('PSD [dBm/Hz]');
drawnow;

%% define fit function and find initial guess
% Here the simulation based function and the intial guess are defined. In
% principle the simulating part "makeDSH_PSD" can use 6 (or even more)
% input parameters for the optimization. However it is strongly recommended
% to keep this number as low as possible. In this example, the center
% frequency and the delay are not fitted, but fixed, since they can be
% predetermined quite well. The remaining vector fp for the fit parameters
% parametrizes the signal strength and the noise sources, x are the
% frequencies. The noises are normalized upon the input and some additional
% inputs were fixed here, in order to feed through adjustments that might
% be made to this script. 
fitFun = @(fp,x) makeDSH_PSD(x,[beatFreq,addedpath,fp(1), ...
    noiseAmplitudeNormalizationFactor*fp(2:4)],'numSamp',numSamp, ...
    'ts',1/(signalSampling*beatFreq));

% Next, the initial guess is made. You should adjust this each and every
% time you use a new dataset to ensure a fast convergence. After adjusting
% the parameters continue in this section to generate a plot in the figure
% window of the already plotted data for comparison. If it does not seem
% close to the data - adjust the fp0 vector and do it again until you are
% satisfied with the result. 
% The parametrization of the vector is as follows:
% fp0 = [signal amplitude, gaussian noise amplitude, flicker noise
% amplitude, random walk noise amplitude];
% 
% here are 'good' guesses for the respective examples (uncomment one or the
% other): 
% ----------------------------------------------------
% example file 14.1.csv
fp0 = [1,15,60,7e5];
% ----------------------------------------------------
% example file 23.1.csv
% fp0 = [1.3,12,60,7e5];
% ----------------------------------------------------
%
% PROCEDURE TO FIND AN ALMOST FITTING INITIAL GUESS 
% (exemplarily for 14.1.csv here ): 
% 1. set fp0's third and fourth entry to 0 
%       >> fp0 = [1,1,0,0];
% 2. adjust first and second entry until the height of the center peak and
% the height of the modulation in the outer area fits 
%       >> fp0 = [1,15,0,0];
% 3. increase the third value until the inner modulation parts start to
% overlap better
%       >> fp0 = [1,15,60,0];
% 4. increase the fourth value until the shape of the central peak fits
%       >> fp0 = [1,15,60,7e5];
% Any optimization by hand is very likely to be much faster than in the
% following automated optimization.

% here the initial guess is plotted - check if it fits already to the data,
% if not, DO NOT PROCEED, but adjust your initial guess fp0 until you have
% a reasonable start point for the optimization 
% -----------------------------------------------------------------------
% FOR CONVERGENCE IT IS VITAL TO HAVE AN ALMOST PERFECTLY FITTING INITIAL GUESS
% -----------------------------------------------------------------------
% The reason for that is the noise on the simulated data that will lead to
% very small (and even varying) gradients for most of the fit parameters
% upon the small adjustments that the fit routine does in order to find the
% direction to alter the parameters to. Being in an optimization-wise
% rather flat area away from optimum is therefore killing the convergence
% of the routine.
% Normal fluctuations of the value of the  summed squared deviation of the
% logarithmic data can be up to 10% for different evaluations of the
% simulation due to the simulations' noise (so optimization steps need to
% be better than that). 
plot(freqs,10*log10(fitFun(fp0,freqs)),'--','Color',[0.5 0.5 0.5]); 
drawnow;

%% fit the simulation to the data using simplex search (Nelder-Mead)
% To fit the simulation to the data, a Nelder-Mead based simplex search
% method is employed. This minimizes the amount of required function
% evaluations (since the simulations take some time) and still ensures
% relatively good convergence to a propper fit. 
% Since it is important to fit especially also the behavior away from the
% main peak the fit is done in logspace.

% Next the options for the fminsearch routine are set. The second line
% there enables the display of the optimization iterations - if this is not
% desired, leave this last options. 
% 
% REMARK: It does not make sense to arbitrarily decrease the required
% accuracy, since the simulation is also noisy and will not fit the data
% increasingly well, if just the target accuracy is increased. To reach
% better values the averaging in the simulation (basically the 'numSamp'
% variable) would need to be increased, which on the other hand linearly
% increases the simulation time, which is the time of just a single
% function evaluation in the optimization. Increased accuracy comes
% therefore at the cost of much larger fit times. 
% If the simulation function is changed (say to a different 'numSamp'
% value) the initial guess will most probably have to be readjusted!
options = optimset('MaxFunEvals',1e2,'TolX',1e-2,'TolFun',1e-2, ...
    'Display','iter');

% In the following the optimization is performed and the optimized
% parameters are put into fp1.
fp1 = fminsearch(@(fp) sum(abs(10*log10(fitFun(fp,freqs))-powerdbm).^2),...
    fp0,options);

% The optimized set is now plotted into again the same window for
% comparison. 
plot(freqs,10*log10(fitFun(fp1,freqs))); hold off
legend({'data','initial guess','NM-optimized'})

%% extract the different noise and relate it to linewidths
% In this part the extracted noise amplitudes for the different noise
% sources are used to generate PSDs and fit the corresponding power law to
% them to extract the needed noise coefficients. 
% The corresponding linewidths are extracted from that. 

% extract the noise coefficients
[kw,kf,kr] = extractLinewidths(noiseAmplitudeNormalizationFactor*fp1(2:4), ...
    beatFreq,'numSamp',numSamp,'ts',1/(signalSampling*beatFreq));

% give the estimated linewidths 
disp(['Estimated lorentian linewidth: ' num2str(kw*pi/1e3) ' kHz'])
tau0 = addedpath/c;
disp(['Estimated gaussian FWHM: ' num2str((sqrt(8*pi^2*kf*log(2)*...
    (4.3+log(17.2*pi^2*kf*tau0^2.1)))/pi)/1e3) ' kHz']) 

% check validity for the approximation of the gaussian linewidth part
if 4*pi^2*kf*tau0^2<=1
    disp('Gaussian FWHM is not a good approximation.')
end
