%%
PhysicsConstants;
InitialLongitudinalTrapFrequency = 56; %in kHz
deltaE = PlanckConstant * InitialLongitudinalTrapFrequency;
FractionOfInitialPotential = 0.01;

GroundStatePopulations = 0.3:0.1:0.9;
InitialTrapDepth = LatticeProperties.estimateTrapDepthFromHeatingSidebandFreqDT1DT3(InitialLongitudinalTrapFrequency);
TrapDepths = linspace(1,FractionOfInitialPotential,1000) .* InitialTrapDepth;
results = LatticeProperties.estimateTrapFreq33FromTrapDepth(0,0,0,0,TrapDepths);
LongitudinalTrapFrequencies = results(1:length(TrapDepths));
TransverseTrapFrequencies = results(length(TrapDepths)+1:end-3);
%LongitudinalTrapFrequencies_Estimated = sqrt(TrapDepths./InitialTrapDepth)*InitialLongitudinalTrapFrequency ;  
%plot(TrapDepths, LongitudinalTrapFrequencies ./ LongitudinalTrapFrequencies_Estimated)
DeltaPClassical = zeros(length(GroundStatePopulations), length(LongitudinalTrapFrequencies));
DeltaPQuantumMechanical = zeros(length(LongitudinalTrapFrequencies), 1);
ltext = {};
for ii = 1:length(GroundStatePopulations)
    InitialTemperature = -(7.24297e22 * deltaE) / log(1 - GroundStatePopulations(ii));
    InitialEntropy = (PlanckConstant*InitialLongitudinalTrapFrequency) / (BoltzmannConstant*InitialTemperature);
    for  jj = 1:length(LongitudinalTrapFrequencies)
        Temperature = (PlanckConstant/(BoltzmannConstant*InitialEntropy)) * LongitudinalTrapFrequencies(jj);
        DeltaPClassical(ii, jj) = Cs133Mass * sqrt(BoltzmannConstant*Temperature/Cs133Mass);
    end
    ltext{end+1} = ['\Delta v_{Classical}--> ' num2str(100*GroundStatePopulations(ii)) '% GS Population'];
end
clear ii jj
for  ii = 1:length(LongitudinalTrapFrequencies)
    DeltaPQuantumMechanical(ii) = sqrt(0.5*Cs133Mass*PlanckConstant*LongitudinalTrapFrequencies(ii));
end
ltext{end+1} = '\Delta v_{QM}';
clear ii
figure(10) 
clf
plot(TrapDepths, DeltaPClassical.*1e6./Cs133Mass)
hold on
plot(TrapDepths, DeltaPQuantumMechanical.*1e6./Cs133Mass, '--')
ylim([0 max(DeltaPClassical(:).*1e6./Cs133Mass)])
sgtitle(['Velocity Spread (1D) v Release Depth (Upto ' num2str(100*FractionOfInitialPotential) '% of initial depth)']);
xlabel('Trap Depth (in kHz)')
ylabel('Uncertainty in Velocity (\mum/s)')
legend(ltext, 'Location', 'NorthWest')