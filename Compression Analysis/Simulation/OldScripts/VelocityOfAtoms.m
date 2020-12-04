%
% All quantities are in SI units
%   Distance --- in metres
%   Time-------- in seconds
%
% Mean velocity in 1D is the average absolute value, since it is 
% zero otherwise due to the symmetric 1-D M-B distribution.
%
% Distance travelled is calculated with Most Probable Velocity
%
%% Initialization
PhysicsConstants;
Lattice = 'HDT';
TrapFrequencyToUse = 'L'; %'L' for Longitudinal or 'T' for Transverse
GroundStatePopulation = 0.8;
n = 1; % Dimensionality
TimeOfFlight = 1e-6;
%%
switch Lattice
    case 'HDT'
        TransverseTrappingFrequencyinHz   = 32e+03;
        LongitudinalTrappingFrequencyinHz = 56e+03;
    case 'VDT'
        Wavelength = 1064e-9;                          % VDT wavelength
        Frequency  = 2*pi*SpeedOfLight/Wavelength;
        detuningD1 = 2*pi*SpeedOfLight*(1/CsD1lambda-1/Wavelength);
        detuningD2 = 2*pi*SpeedOfLight*(1/CsD2lambda-1/Wavelength);
        k  = 2*pi./Wavelength;
        w0 = 50e-6;                                    % Beam waist
        P  = 70e-3;                                    % Beam power
        fOscD1 = 0.344;                                % D1 Absorption oscillator strength
        fOscD2 = 0.714;                                % Absorption oscillator strength
        potentialContributionD1 = (fOscD1*CsD1Gamma/(2*pi*SpeedOfLight/CsD1lambda)^3)*(1/detuningD1+(1/(detuningD1+2*Frequency)));
        potentialContributionD2 = (fOscD2*CsD2Gamma/(2*pi*SpeedOfLight/CsD2lambda)^3)*(1/detuningD2+(1/(detuningD2+2*Frequency)));
        I0 = 2.*P./(pi.*(w0).^2);  % Single peak beam intensity, factor of 2 because two counter propagating beams
        Imax = 2*I0; % Multiply another factor of 2 because of interference
        U0 = -(3*pi*SpeedOfLight^2*Imax/2)*(potentialContributionD1+potentialContributionD2);
        TransverseTrappingFrequencyinHz = (1/(2*pi)) *  sqrt(4*abs(U0)/(Cs133Mass*w0^2));
        LongitudinalTrappingFrequencyinHz = sqrt(2*abs(U0)/(Cs133Mass*Wavelength^2));
end
%% Temperature of Atoms given Trap Frequency and GS population
switch TrapFrequencyToUse
    case 'T'
        TrapFrequency = TransverseTrappingFrequencyinHz;
    case 'L'
        TrapFrequency = LongitudinalTrappingFrequencyinHz;
end
deltaE = PlanckConstant * TrapFrequency;
Temperature = -(7.24297e22 * deltaE) / log(1 - GroundStatePopulation);

%% Calculation of velocities
MostProbableVelocity = sqrt(((n-1)*BoltzmannConstant*Temperature)/Cs133Mass)
MeanVelocity = sqrt((2*BoltzmannConstant*Temperature)/Cs133Mass) * gamma((n+1)/2)/gamma(n/2)
RMSVelocity = sqrt((n*BoltzmannConstant*Temperature)/Cs133Mass)
DistanceTravelled = MostProbableVelocity * TimeOfFlight

    