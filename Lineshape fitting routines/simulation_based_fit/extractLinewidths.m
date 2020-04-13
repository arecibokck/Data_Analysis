function [kw,kf,kr] = extractLinewidths(noiseAmplitudes,beatFreq,varargin)
% This function can be used to extract the different noise spectra
% parameters (see e.g. https://doi.org/10.1364/AO.58.003555 and
% https://doi.org/10.1109/JLT.2008.925046) from a fit of the noise
% amplitudes from a double sided heterodyne measurement. For more see
% simuFitScript.m
%
% call:
% [kw,kf,kr] = extractLinewidths(noiseAmplitudes,beatFreq)
%
% ... to modify internal variables an additional variable input in the form
% of extractLinewidths(...,'variableName',value) can be made.
%
% ---------- Output ----------
% kw                    single sided noise coefficient for gaussian (white)
%                       noise [Hz^2/Hz] (1/f^2 for S_{\phi\phi})
% kf                    single sided noise coefficient of the flicker noise
% kr                    single sided noise coefficient for the random walk
% 
% ---------- Input -----------
% noiseAmplitudes       noise amplitudes as for example used in makeDSH_PSD
%                       in conjunction with frequencyNoise2.m to generate
%                       phase noise
% beatFreq              frequency constant used for the sampling (sampling
%                       time units are: 1/(5*beatFreq)) - same as in
%                       makeDSH_PSD 
%
% varargin: some internal variables can be modified via additional input
% pairs of the form 'variable name',value - in the following, the default
% values are shown: 
% 'numSamp',1e7         % this gives the length of the employed time
%                       % vector, use the same as in the code used to find
%                       % the noise amplitudes
% 'ts', ...             % sampling time, default value calculated as
%                       % 1/(5*beatFreq)) - use same as in the code for
%                       % finding noise amplitudes
% 'cutOff',1            % 1/f^n is not defined for f=0 so introduce lower
%                       % cut-off for the considered frequencies

%% initial parameters and potential overwriting
numSamp = 1e7;
ts = 1/(5*beatFreq);
cutOff = 1; % 25/(ts*numSamp);

for Index = 1:numel(varargin)/2
    eval([varargin{2*Index-1} '= varargin{2*Index};']);
end

%% initialize figure
figure; hold on;
set(gcf, 'WindowStyle', 'Docked');
noiseCell = {'kw','kf','kr'};

%% fit each noise separately and determine the respective k_noise
for Index = 1:numel(noiseAmplitudes)
    
    % current noise amplitudes
    currentNoise = zeros(1,3);
    currentNoise(Index) = noiseAmplitudes(Index);
    
    % generate noise sample and S_{\phi\phi}
    [timeVector,phaseNoise] = frequencyNoise2(ts,numSamp,currentNoise);
    [freq0,phaseNoiseSpecDens] = ownWelch(timeVector,phaseNoise);
    
    % get rid of DC stuff (otherwise error for f = 0)
    phaseNoiseSpecDens(freq0<cutOff) = [];
    freq0(freq0<cutOff) = [];
    
    % go to log-data to improve the fit
    logDataX = log10(freq0);
    logDataY = log10(phaseNoiseSpecDens);
    
    % get a good guess 
    midIndex = round(numel(freq0)/2);
    fp0 = [(phaseNoiseSpecDens(midIndex)*freq0(midIndex)^(Index+1)),phaseNoiseSpecDens(end)];
    
    % current function to fit and fit
    phaseFitFun = @(a,b,x) log10(a*fp0(1)*x.^(-(Index+1))+fp0(2)*b);
    fitt = fittype(phaseFitFun);
    fitParam = fit(freq0',logDataY',fitt,'StartPoint',[1 1],'Lower',[0,0]);

    % extract fitted parameters
    fp = coeffvalues(fitParam).*fp0;
    fP = coeffvalues(fitParam);
    
    % plot them
    plot(logDataX,logDataY);
    plot(logDataX,phaseFitFun(fP(1),fP(2),freq0));
    xlabel('log10(\nu)'); ylabel('log10(S_{\phi\phi})');
    drawnow;
    
    % evaluate the current noise parameter
    eval([noiseCell{Index} '=' num2str(fp(1)) ';']);
end

%% last labels of the fit
title('phase noise fits');
legend({'gaussian',['k_w = ' num2str(kw)],'flicker',['k_f = ' num2str(kf)], ...
    'random',['k_r = ' num2str(kr)],});

end