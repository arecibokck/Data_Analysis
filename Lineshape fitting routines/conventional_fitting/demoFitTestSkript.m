clear all
close all

% This is a test script to read in and fit data with the conventional Voigt
% fitting for an electrical spectrum analyzer (ESA) measurement of a
% delayed self heterodyning setup for the linewidth measurement of a laser.

%% read data
fileName = 'demoSpectrum_14.csv';
data = csvread(fileName);
freqs = data(:,1);
powerdbm = data(:,2);
delF = mean(freqs);

%% define fit function and guess
fitFun = @(fp,x) generatePSD(x-fp(6),fp); % documentation see there

% guess for: [linewidth/(2*pi), delay time, power scaling shot noise, power
% scaling Voigt, background noise, frequency displacement] 
fp0 = [20e3/(2*pi), 0.021e-3, 1e6, 1e5, 1e-11, delF]; 

%% do the fit
% set least square fit options
options = optimoptions('lsqcurvefit','StepTolerance',1e-8,...
    'MaxFunctionEvaluations',5e3,'MaxIterations',1e3);
% do least square fit
FP = lsqcurvefit(fitFun,fp0,freqs,powerdbm,[],[],options);

% set Nelder-Mead options
options = optimset('TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',1e4,'MaxIter',1e4);
% do fminsearch
FP2 = fminsearch(@(fp) sum(abs(fitFun(fp,freqs)-powerdbm).^2),...
    fp0,options);

% give estimated linewidths
disp('The estimated linewidths (with lsq and NM method)'); 
disp(['are ' num2str(2*pi*[FP(1) FP2(1)]*1e-3) ' kHz.']);

%% plotting
% plot lsq fit
plot(freqs,fitFun(FP,freqs)); hold on;
% plot NM-fit
plot(freqs,fitFun(FP2,freqs));
% plot data
plot(freqs,powerdbm,'.');  
% plot initial guess
plot(freqs,fitFun(fp0,freqs),'-.','Color',[0.8 0.8 0.8]); 

% proper labeling:
xlabel('Beat frequency [Hz]');
ylabel('PSD [dBm/Hz]');
legend({'least square fit','Nelder-Mead "fit"', 'ESA data', 'initial guess'});

% rearrange ...
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));

% save stuff
saveas(gcf,['fit_' fileName(1:end-3) 'fig']);
saveas(gcf,['fit_' fileName(1:end-3) 'png']);


