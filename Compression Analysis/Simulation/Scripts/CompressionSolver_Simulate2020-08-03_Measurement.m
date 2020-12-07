%% - Create solver object with specified options
clc
%%
options={};
%options=[options,{'Parameter','Value'}];
options=[options,{'TimeStep',                         10e-06}];   
options=[options,{'FinalTime',                       2.5e-03}];   
options=[options,{'VDTPower',                         500e-3}];
options=[options,{'HDTPower',                          12e-3}];   
options=[options,{'PotentialType',                'gaussian'}];   
options=[options,{'CompressionDirection',       'horizontal'}]; 
options=[options,{'PositionDistributionLimits', [-1 1]*60E-6}];   
options=[options,{'VelocityDistributionLimits', [-1 1]*8E-3}]; 
options=[options,{'PositionSampling', 0.1E-6}]; 
options=[options,{'VelocitySampling', 0.05E-3}]; 

Solver  = CompressionSolver(options{:});

%% - plot potential along with its harmonic approximation and associated forces as a validation
Solver.plotAllPotentials;

%% - specify initial conditions for simulation
OptionsStruct = struct;
OptionsStruct.GroundStatePopulation                            = 0.9;
OptionsStruct.FinalTrapDepthInUnitsOfRecoilEnergy                = 5;
OptionsStruct.TypeOfPositionDistribution = 'mirroredFlatTopGaussian';
OptionsStruct.MeanForPositionDistribution                    = -0e-6;
OptionsStruct.SDForPositionDistribution                      = 15e-6;
OptionsStruct.TypeOfVelocityDistribution       = 'maxwell-boltzmann';
OptionsStruct.NumberOfAtoms                                  = 1000;
OptionsStruct.FlatWidth                                      = 19e-6;

options = CompressionHelper.convertstruct2cell(OptionsStruct);

[initialPositions, initialVelocities, PositionDistribution, VelocityDistribution] = Solver.setInitialConditions(options{:});

Solver.plotInitialConditions('NumberOfBins',50);

%% - Plot initial conditions as sampled along with the distribution they were sampled from
NumberOfBins = 100;
Plotting.plotSampling(Solver, PositionDistribution, VelocityDistribution, NumberOfBins);

%% - run simulation and Time evolution in Position Space
poolobj = gcp('nocreate'); % Check if pool is open
if isempty(poolobj)
    parpool;
end
tic
Solver.runClassicalSimulation()
toc
% - 
TimeForFullCompression = Plotting.plotTrajectories(Solver);

% - Extracting quarter periods from the trajectories
Plotting.plotQuarterPeriods(Solver, TimeForFullCompression)

% - Time evolution of RMS Spread
Plotting.plotSpreadEvolution(Solver)

%% - Time evolution in Phase Space
NumberOfBins = 100;
Plotting.plotPhaseSpaceEvolution(Solver, NumberOfBins)

%% - Analysis of time evolution in Phase Space

%% - Time evolution of RMS Spread for Different Temperatures
Spreads = {};
GroundStatePopulations = 0.3:0.2:0.9;
OptionsStruct = Solver.InitialDistributionParameters;
for ii = 1:numel(GroundStatePopulations)
    OptionsStruct.GroundStatePopulation = GroundStatePopulations(ii);
    options = CompressionHelper.convertstruct2cell(OptionsStruct);
    [initialPositions, initialVelocities] = Solver.setInitialConditions(options{:});
    Solver.runClassicalSimulation()
    positions = Solver.simulationResults(:,:,1);
    RMSSpread = zeros(size(positions, 1),1);
    for Index = 1:size(positions, 1)
        RMSSpread(Index) = rms(positions(Index,:));
    end
    Spreads{end+1} = RMSSpread;
end

Plotting.plotSpreadEvolutionForDifferentTemps(Solver, GroundStatePopulations, Spreads)

%% - First minimum spread for different initial spreads of distributions and temperatures
tic

MinimumSpreads = {};
PeriodsAtMinimumSpread = {};
GroundStatePopulations = [0.3 0.7 0.9];
OptionsStruct = Solver.InitialDistributionParameters;
for ii = 1:numel(GroundStatePopulations)
    InitialSpreads = 1:1:25;
    FirstMinimumSpreads = zeros(size(InitialSpreads, 1),1);
    CorrespondingPeriod = zeros(size(InitialSpreads, 1),1);
    OptionsStruct.GroundStatePopulation = GroundStatePopulations(ii);
    
    for ii = 1:length(InitialSpreads)
        OptionsStruct.SDForPositionDistribution = InitialSpreads(ii) * 1e-6;
        options = CompressionHelper.convertstruct2cell(OptionsStruct);
        [initialPositions, initialVelocities] = Solver.setInitialConditions(options{:});
        Solver.runClassicalSimulation()
        positions = Solver.simulationResults(:,:,1);
        RMSSpread = zeros(size(positions, 1),1);
        for Index = 1:size(positions, 1)
            RMSSpread(Index) = rms(positions(Index,:));
        end
        [FirstMinimumSpreads(ii), idx] = min(RMSSpread);
        CorrespondingPeriod(ii) = Solver.timeSpan(idx);
    end
    MinimumSpreads{end+1} = FirstMinimumSpreads;
    PeriodsAtMinimumSpread{end+1} = CorrespondingPeriod;
end

Plotting.plotMinimumvInitialSpreads(Solver, InitialSpreads, MinimumSpreads, GroundStatePopulations)

toc
%% - Time evolution of RMS Spread for Different Release Depths

%% - Cycling of Compression


