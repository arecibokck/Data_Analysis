function reinitializeCompressionSolver(this,positionSampling,velocitySampling)
    %% PHYSICAL CONSTANTS
    pc = this.physicalConstants;
    %% TRAP FREQUENCIES & DEPTHS
    % trap frequencies and potential depth of VDT
    Frequency                     = 2*pi*pc.SpeedOfLight/this.VDTWavelength;
    detuningD1                    = 2*pi*pc.SpeedOfLight*(1/pc.CsD1lambda-1/this.VDTWavelength);
    detuningD2                    = 2*pi*pc.SpeedOfLight*(1/pc.CsD2lambda-1/this.VDTWavelength);
    fOscD1                        = 0.344;                                        % D1 Absorption oscillator strength
    fOscD2                        = 0.714;                                        % Absorption oscillator strength
    potentialContributionD1       = (fOscD1*pc.CsD1Gamma/(2*pi*pc.SpeedOfLight/pc.CsD1lambda)^3)*(1/detuningD1+(1/(detuningD1+2*Frequency)));
    potentialContributionD2       = (fOscD2*pc.CsD2Gamma/(2*pi*pc.SpeedOfLight/pc.CsD2lambda)^3)*(1/detuningD2+(1/(detuningD2+2*Frequency)));
    I0                            = 2.*this.VDTPower./(pi.*(this.VDTBeamWaist).^2);  % Single peak beam intensity, factor of 2 because two counter propagating beams
    Imax                          = 2*I0; % Multiply another factor of 2 because of interference
    this.potentialDepthVDT        = -(3*pi*pc.SpeedOfLight^2*Imax/2)*(potentialContributionD1+potentialContributionD2); % (eq. 1.22a Thesis alt)
    trapFreqVDTLongitudinalinHz   = sqrt(2*abs(this.potentialDepthVDT)/(pc.Cs133Mass*this.VDTWavelength^2)); % (eq. 1.45 Thesis W.Alt)
    trapFreqVDTTransverseinHz     = (1/(2*pi)) *  sqrt(4*abs(this.potentialDepthVDT)/(pc.Cs133Mass*this.VDTBeamWaist^2)); % (eq. 1.46 Thesis W.Alt)
    this.trapFreqVDTLongitudinal  = (2*pi) * trapFreqVDTLongitudinalinHz;
    this.trapFreqVDTTransverse    = (2*pi) * trapFreqVDTTransverseinHz;
    % trap frequencies and potential depth of HDT
    trapFreqHDTyinHz              = sqrt((this.HDTPower/12e-3)) * 56e+03; % Known trap frequency from heating sideband
    this.potentialDepthHDT        = - 1e3 * pc.PlanckConstant * LatticeProperties.estimateTrapDepthFromHeatingSidebandFreqDT1DT3(trapFreqHDTyinHz*1e-3);
    result                        = LatticeProperties.estimateTrapFreq33FromTrapDepth(0,0,0,0,(abs(this.potentialDepthHDT)/pc.PlanckConstant)*1e-3);
    trapFreqHDTxinHz              = result(2) * 1e+03;
    trapFreqHDTzinHz              = sqrt((this.HDTPower/12e-3)) * 0.534e+03;
    this.trapFreqHDTx             = (2*pi) * trapFreqHDTxinHz;
    this.trapFreqHDTy             = (2*pi) * trapFreqHDTyinHz;
    this.trapFreqHDTz             = (2*pi) * trapFreqHDTzinHz;
    
    %% SIMULATION PARAMETERS
    this.positionGrid   = this.positionLimits.Lower :positionSampling: this.positionLimits.Upper;
    this.velocityGrid   = this.velocityLimits.Lower :velocitySampling: this.velocityLimits.Upper;
end