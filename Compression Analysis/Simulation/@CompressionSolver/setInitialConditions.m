function [initialPositions, initialVelocities, varargout] = setInitialConditions(this, varargin)
    
    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'GroundStatePopulation', 0.9,...
            @(x) assert(isnumeric(x) && isscalar(x) && (x > 0) && (x < 1)));  
    addParameter(p, 'FinalTrapDepthInUnitsOfRecoilEnergy', 5,...
            @(x) assert(isnumeric(x) && isscalar(x) && (x >= 0)));  
    addParameter(p, 'TypeOfPositionDistribution', 'gaussian',...
            @(x) assert(any(strcmpi(x,{'gaussian','higherOrderGaussian', 'mirroredFlatTopGaussian'})), 'Incorrect entry! Check and try again!'));  
    addParameter(p, 'MeanForPositionDistribution', 0, ...
            @(x) assert(isnumeric(x) && isscalar(x)));  
    addParameter(p, 'SDForPositionDistribution',   5e-06, ...
            @(x) assert(isnumeric(x) && isscalar(x)));  
    addParameter(p, 'TypeOfVelocityDistribution', 'maxwell-boltzmann', ...
            @(x) assert(any(strcmpi(x,{'zero', 'uniform', 'gaussian', 'maxwell-boltzmann'})), 'Incorrect entry! Check and try again!'));
    addParameter(p, 'MeanForVelocityDistribution', 0,...
            @(x) assert(isnumeric(x) && isscalar(x)));  
    addParameter(p, 'SDForVelocityDistribution',   5e-06,...
            @(x) assert(isnumeric(x) && isscalar(x)));  
    addParameter(p, 'NumberOfAtoms',   10000,...
            @(x) assert(isnumeric(x) && isscalar(x) && (x > 0)));
    addParameter(p,'FlatWidth',10e-6,...
            @(x) x>=0 && isscalar(x));
    p.parse(varargin{:});
    
    this.numberOfAtoms      = p.Results.NumberOfAtoms;
        
    GroundStatePopulation = p.Results.GroundStatePopulation;
    FinalTrapDepthInUnitsOfRecoilEnergy = p.Results.FinalTrapDepthInUnitsOfRecoilEnergy;
    TypeOfPositionDistribution = p.Results.TypeOfPositionDistribution;
    TypeOfVelocityDistribution = p.Results.TypeOfVelocityDistribution;
    flatWidth = p.Results.FlatWidth;
    InitialTrapDepthInUnitsOfRecoilEnergy = abs(this.potentialDepth)/this.recoilEnergy;
    
    
    pc = this.physicalConstants;
    %PowerAtTargetReleaseDepth = (((2 * FinalTrapDepthInUnitsOfRecoilEnergy * this.recoilEnergy)/(3*pi * pc.SpeedOfLight^2 * (potentialContributionD1+potentialContributionD2))) * (pi.*(this.beamWaist).^2))/2;
    FractionOfInitialPotential = FinalTrapDepthInUnitsOfRecoilEnergy / InitialTrapDepthInUnitsOfRecoilEnergy;
    deltaE = pc.PlanckConstantReduced * this.relevantTFForTemperature;
    initialTemperatureBeforeAdiabaticRampDown = -(7.24297e22 * deltaE) / log(1 - GroundStatePopulation);  % See sec. 2 of note for derivation:  https://www.evernote.com/shard/s124/nl/13819157/a7d81880-899b-426f-82e6-f8d10e6c094b
    initialTemperature = sqrt(FractionOfInitialPotential) * initialTemperatureBeforeAdiabaticRampDown; % See sec. 2 of note for derivation:  https://www.evernote.com/shard/s124/nl/13819157/a7d81880-899b-426f-82e6-f8d10e6c094b
    
    n = this.numberOfAtoms;
    
    % - sampling the position distribution
    positions = this.positionGrid;
    mean = p.Results.MeanForPositionDistribution;
    sd = p.Results.SDForPositionDistribution;
    switch TypeOfPositionDistribution
        case 'gaussian'
            PositionDistribution  = CompressionHelper.GaussianDistribution(mean, (sd)^2, positions); % Mean, Variance, Position Vector
            initialPositions = CompressionHelper.drawSamplesFromDistribution(n, positions, PositionDistribution); 
        case 'higherOrderGaussian'
            PositionDistribution  = CompressionHelper.HigherOrderGaussianDistribution(mean, (sd)^2, 2, positions); % Mean, Variance, Position Vector
            initialPositions = CompressionHelper.drawSamplesFromDistribution(n, positions, PositionDistribution); 
        case 'mirroredFlatTopGaussian' 
            PositionDistribution  = CompressionHelper.MirroredFlatTopGaussianDistribution(mean, flatWidth, sd, positions); % Mean, Variance, Position Vector
            initialPositions = CompressionHelper.drawSamplesFromDistribution(n, positions, PositionDistribution); 
    end
    
    % - sampling the velocity distribution
    velocities = this.velocityGrid;
    mean = p.Results.MeanForVelocityDistribution;
    sd   = p.Results.SDForVelocityDistribution;
    switch TypeOfVelocityDistribution
        case 'zero'
            initialVelocities  = 0;
        case 'uniform'
            VelocityDistribution     = CompressionHelper.UniformDistribution(mean, sd, velocities);
            initialVelocities = CompressionHelper.drawSamplesFromDistribution(n, velocities, VelocityDistribution);
        case 'gaussian' 
            VelocityDistribution     = CompressionHelper.GaussianDistribution(mean, sd, velocities);
            initialVelocities = CompressionHelper.drawSamplesFromDistribution(n, velocities, VelocityDistribution);
        case 'maxwell-boltzmann' 
            VelocityDistribution     = CompressionHelper.MaxwellBoltzmannDistribution(initialTemperature, velocities);
            initialVelocities = CompressionHelper.drawSamplesFromDistribution(n, velocities, VelocityDistribution);
    end
    signflips = (-1) .* (rand(length(initialVelocities),1) > 0.5);
    signflips(signflips == 0) = 1;
    initialVelocities = initialVelocities .* signflips;
    %
    this.initialPositions = initialPositions;
    this.initialVelocities = initialVelocities;
    
    % - output handling
    if nargout >0
    varargout = {PositionDistribution, VelocityDistribution};
    end
    
    % - store in struct
    this.InitialDistributionParameters = struct;
    this.InitialDistributionParameters.GroundStatePopulation = GroundStatePopulation;
    this.InitialDistributionParameters.FinalTrapDepthInUnitsOfRecoilEnergy = FinalTrapDepthInUnitsOfRecoilEnergy;
    this.InitialDistributionParameters.TypeOfPositionDistribution = TypeOfPositionDistribution;
    this.InitialDistributionParameters.TypeOfVelocityDistribution = TypeOfVelocityDistribution;
    this.InitialDistributionParameters.FractionOfInitialPotential = FractionOfInitialPotential;
    this.InitialDistributionParameters.initialTemperature = initialTemperature;
    this.InitialDistributionParameters.flatWidth = flatWidth;
end