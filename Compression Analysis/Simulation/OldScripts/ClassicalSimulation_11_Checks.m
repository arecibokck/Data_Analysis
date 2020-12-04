clearvars
clc
poolobj = gcp('nocreate'); % Check if pool is open
if isempty(poolobj)
    parpool;
end
DebugMode = false;
%% initializing the VDT potential
PhysicsConstants;
Wavelength = 1064e-9;                               % VDT wavelength
Frequency  = 2*pi*SpeedOfLight/Wavelength;
detuningD1 = 2*pi*SpeedOfLight*(1/CsD1lambda-1/Wavelength);
detuningD2 = 2*pi*SpeedOfLight*(1/CsD2lambda-1/Wavelength);
% Trap parameters
Trap.k = 2*pi./Wavelength;
Trap.w0   = 50e-6;                                  % Beam waist
Trap.z_R  = pi.*(Trap.w0.^2)/Wavelength;            % Rayleigh length
Trap.w = @(z) Trap.w0.*sqrt(1+((z./Trap.z_R).^2));  % Axial waist
Trap.P = 500e-3;                                    % Beam power
fOscD1 = 0.344;                                     % D1 Absorption oscillator strength
fOscD2 = 0.714;                                     % Absorption oscillator strength
potentialContributionD1 = (fOscD1*CsD1Gamma/(2*pi*SpeedOfLight/CsD1lambda)^3)*(1/detuningD1+(1/(detuningD1+2*Frequency)));
potentialContributionD2 = (fOscD2*CsD2Gamma/(2*pi*SpeedOfLight/CsD2lambda)^3)*(1/detuningD2+(1/(detuningD2+2*Frequency)));
I0 = 2.*Trap.P./(pi.*(Trap.w0).^2);  % Single peak beam intensity, factor of 2 because two counter propagating beams
Imax = 2*I0; % Multiply another factor of 2 because of interference
Trap.U0 = -(3*pi*SpeedOfLight^2*Imax/2)*(potentialContributionD1+potentialContributionD2);
Trap.U0InTemperature = Trap.U0/BoltzmannConstant;
Trap.U0InFreq = Trap.U0/PlanckConstant;
%% Temperature induced broadening of velocity (momentum) distribution post adiabatic release by ramping down the HDT for horiontal compression in the VDT
LongitudinalTrapFrequency = 56; %in kHz
GroundStatePopulation = 0.3;
FractionsOfInitialPotential = linspace(0.01, 1, 100);
initialTemperature = zeros(length(FractionsOfInitialPotential), 1);
for Index = 1:length(FractionsOfInitialPotential)
    deltaE = PlanckConstantReduced * 2 * pi * LongitudinalTrapFrequency;
    initialTemperatureBeforeAdiabaticRampDown = -(7.24297e22 * deltaE) / log(1 - GroundStatePopulation);
    initialTemperature(Index) = sqrt(FractionsOfInitialPotential(Index)) * initialTemperatureBeforeAdiabaticRampDown;
end

GroundStatePopulations = linspace(0.1, 0.9, 100);
initialTemperatureBeforeAdiabaticRampDown = zeros(length(GroundStatePopulations), 1);
for Index = 1:length(GroundStatePopulations)
    deltaE = PlanckConstantReduced * 2 * pi * LongitudinalTrapFrequency;
    initialTemperatureBeforeAdiabaticRampDown(Index) = -(7.24297e22 * deltaE) / log(1 - GroundStatePopulations(Index));
end

plotChecks(GroundStatePopulations, initialTemperatureBeforeAdiabaticRampDown, FractionsOfInitialPotential, initialTemperature)

function plotChecks(GroundStatePopulations, initialTemperatureBeforeAdiabaticRampDown, FractionsOfInitialPotential, initialTemperature)
% - plot Distribution
figure(11)
clf
ProbDFName = {'Temp. v GS-Population','Temp. v Fraction of Initial Potential'};
for kk = 1:4
    subplot(2,2,kk)
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    set(gcf, 'defaultAxesColorOrder', [[0 0 0];[0 0 0]]);
    if kk == 1
        plot(GroundStatePopulations, initialTemperatureBeforeAdiabaticRampDown*1e9, 'Color', colours{1})
        xlabel('GS Populations','FontSize', 14)
        ylabel('Temperature (nK)','FontSize', 14)
        title(ProbDFName{1},'FontSize', 14)
    elseif kk == 2
        plot(FractionsOfInitialPotential, initialTemperature*1e9, 'Color', colours{2})
        xlabel('Fractions of Initial Potential','FontSize', 14)
        ylabel('Temperature (nK)','FontSize', 14)
        legend('GS Pop = 0.3')
        title(ProbDFName{2},'FontSize', 14)
    elseif kk == 3
        PhysicsConstants;
        NumberOfAtoms = 1000;
        max_vel = 2e-3;
        velocities = 0:(max_vel/(NumberOfAtoms-1)):max_vel;
        LongitudinalTrapFrequency = 56; %in kHz
        FractionOfInitialPotential = 0.05;
        GroundStatePopulations = linspace(0.1, 0.9, 3);
        ltext = {};
        for Index = 1:length(GroundStatePopulations)
            deltaE = PlanckConstantReduced * 2 * pi * LongitudinalTrapFrequency;
            initialTemperatureBeforeAdiabaticRampDown = -(7.24297e22 * deltaE) / log(1 - GroundStatePopulations(Index));
            initialTemperature = sqrt(FractionOfInitialPotential) * initialTemperatureBeforeAdiabaticRampDown;
            MBD = MaxwellBoltzmannDistribution(initialTemperature, velocities);
            MBD     = MBD./sum(MBD);
            MBD = horzcat(flip(MBD), MBD);
            VelocityRange = horzcat(-flip(velocities*1e3), velocities*1e3);
            plot(VelocityRange, MBD, 'Color', colours{Index})
            ltext{end + 1} = ['GS Pop = ' num2str(GroundStatePopulations(Index))];
            hold on
        end
        hold off
        xlabel('Velocities (mm/s)','FontSize', 14)
        ylabel('Probability','FontSize', 14)
        legend(ltext)
    elseif kk == 4
        PhysicsConstants;
        NumberOfAtoms = 1000;
        max_vel = 2e-3;
        velocities = 0:(max_vel/(NumberOfAtoms-1)):max_vel;
        LongitudinalTrapFrequency = 56; %in kHz
        GroundStatePopulation = 0.3;
        FractionsOfInitialPotential = [0.01, 0.05, 0.9];
        ltext = {};
        for Index = 1:length(FractionsOfInitialPotential)
            deltaE = PlanckConstantReduced * 2 * pi * LongitudinalTrapFrequency;
            initialTemperatureBeforeAdiabaticRampDown = -(7.24297e22 * deltaE) / log(1 - GroundStatePopulation);
            initialTemperature = sqrt(FractionsOfInitialPotential(Index)) * initialTemperatureBeforeAdiabaticRampDown;
            MBD = MaxwellBoltzmannDistribution(initialTemperature, velocities);
            MBD     = MBD./sum(MBD);
            MBD = horzcat(flip(MBD), MBD);
            VelocityRange = horzcat(-flip(velocities*1e3), velocities*1e3);
            plot(VelocityRange, MBD, 'Color', colours{Index})
            ltext{end + 1} = ['Fraction = ' num2str(FractionsOfInitialPotential(Index))];
            hold on
        end
        hold off
        xlabel('Velocities (mm/s)','FontSize', 14)
        ylabel('Probability','FontSize', 14)
        lg = legend(ltext);
        title(lg, ['GS Pop = ' num2str(GroundStatePopulation)])
    end 
end
sgtitle('Checks','FontSize', 18)
end
%% helper functions
function ret = UniformDistribution(a, b, x)
    ret = zeros(length(x),1);
    ret(x>a & x<b) = 1/(b - a);
end
function ret = GaussianDistribution(mean, var, x)
    % return the gauss distribution with parameters mean, var.
    % p = gaussdis(mean, var, x);
    ret = 1/sqrt(2*pi*var)*exp(-((x-mean).^2)/(2*var));
end
function ret = MaxwellBoltzmannDistribution(Temperature, Velocity)
    % ret = MaxwellBoltzmannDistribution(Temperature, Velocity)
    % computes the Probability density for a Cs133-atom to have a velocity
    % [m/s] at given Temperature [K]
    
    PhysicsConstants;
    m = Cs133Mass;
    k = BoltzmannConstant;
    %ret = 4*pi*sqrt((m/(2*pi*k*Temperature))^3) ...
    %     * Velocity.^2  .* exp((-m*Velocity.^2)/(2*k*Temperature));
    ret = sqrt((m/(2*pi*k*Temperature))^(1/2)) ...
         .* exp((-m*Velocity.^2)/(2*k*Temperature));
    ret = ret./sum(ret);
end
function ret = drawSamplesFromDistribution(NumberOfSamples, Xvals, ProbabilityDensityFunction)
    % ret = drawSamplesFromDistribution(NumberOfAtoms, velocities, ProbabilityDensityFunction)
    % ret is a N*1 array of randomly sampled values from a distribution defined by P(X) 
    %
    % input: NumberOfSamples   
    %        Xvals   
    %        ProbabilityDensityFunction   
    %
    % output:  ret    N*1 array of randomly sampled values
    %
    % Note:  this implements "Inverse transform sampling"
    
    % - normalize in L1-norm
    ProbabilityDensityFunction     = ProbabilityDensityFunction./sum(ProbabilityDensityFunction); 
    % - compute cumulative distribution
    CumulativeDistributionFunction = cumsum(ProbabilityDensityFunction); 
    % - append a 0 at the beginning of the distribution if the input starts with a non-zero entry (e.g. when sampling a gaussian-distribution)
    if(CumulativeDistributionFunction(1)>0) 
        CumulativeDistributionFunction = [0 CumulativeDistributionFunction];
        dx = Xvals(2) - Xvals(1);
        x0 = Xvals(1) - dx;
        Xvals = [x0 Xvals];
    end 
    % - remove repeted values
    [CumulativeDistributionFunction, index] = unique(CumulativeDistributionFunction); 
    
    % - create uniformly distributed random numbers between 0 and 1
    rnd = rand(NumberOfSamples, 1); 
    % - calculate output
    ret = interp1(CumulativeDistributionFunction, Xvals(index), rnd, 'linear', 0);
end