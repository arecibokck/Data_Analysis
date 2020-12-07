clearvars
clc
poolobj = gcp('nocreate'); % Check if pool is open
if isempty(poolobj)
    parpool;
end
DebugMode = false;
%% initializing the HDT potential
PhysicsConstants;
% Trap parameters
Trap.LongitudinalTrapFrequencyinHz = 56e+03;
Trap.U0 = - 1e3 * PlanckConstant * LatticeProperties.estimateTrapDepthFromHeatingSidebandFreqDT1DT3(Trap.LongitudinalTrapFrequencyinHz*1e-3);
Trap.U0InTemperature = abs(Trap.U0)/BoltzmannConstant;
Trap.U0InFreq = abs(Trap.U0)/PlanckConstant;
%result = LatticeProperties.estimateTrapFreq33FromTrapDepth(0,0,0,0,Trap.U0InFreq*1e-3);
%Trap.TransverseTrapFrequencyinHz = result(2) * 1e3;
Trap.LongitudinalTrapFrequency = (2*pi) *  Trap.LongitudinalTrapFrequencyinHz;
%Trap.TransverseTrapFrequency = (2*pi) *  Trap.TransverseTrapFrequencyinHz;
Trap.w0   = 25e-6;     % Beam waist
%% Temperature induced broadening of velocity (momentum) distribution post adiabatic release by ramping down the HDT for horiontal compression in the VDT
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

% TransverseTrappingFrequencyinHz = (1/(2*pi)) *  sqrt(4*abs(U0)/(Cs133Mass*w0^2));
% TransverseTrappingFrequency = (2*pi) *  TransverseTrappingFrequencyinHz;

LongitudinalTrappingFrequencyinHz = sqrt(2*abs(U0)/(Cs133Mass*Wavelength^2));
LongitudinalTrappingFrequency = (2*pi) *  LongitudinalTrappingFrequencyinHz;

GroundStatePopulation = 0.3;

RecoilEnergy = (PlanckConstantReduced * (2*pi/Wavelength))^2 / (2*Cs133Mass);
InitialTrapDepthInUnitsOfRecoilEnergy = abs(U0)/RecoilEnergy;
FinalTrapDepthInUnitsOfRecoilEnergy = 5; 

PowerAtTargetReleaseDepth = (((2 * FinalTrapDepthInUnitsOfRecoilEnergy * RecoilEnergy)/(3*pi * SpeedOfLight^2 * (potentialContributionD1+potentialContributionD2))) * (pi.*(w0).^2))/2;

FractionOfInitialPotential = FinalTrapDepthInUnitsOfRecoilEnergy / InitialTrapDepthInUnitsOfRecoilEnergy;

deltaE = PlanckConstantReduced * LongitudinalTrappingFrequency;
initialTemperatureBeforeAdiabaticRampDown = -(7.24297e22 * deltaE) / log(1 - GroundStatePopulation);
initialTemperature = sqrt(FractionOfInitialPotential) * initialTemperatureBeforeAdiabaticRampDown;
%% Simulation of trajectory of an atom allowed to oscillate in the trap
NumberOfAtoms = 10000;
max_pos = 20e-6;
positions = linspace(-max_pos,max_pos,NumberOfAtoms);
PositionDistributionTypes = {'InBuiltNormal', 'Gaussian', 'HigherOrderGaussian', 'MirroredFlatTopGaussian'};
ChoiceOfType = 2;

switch ChoiceOfType
    case 1
        initialPositions  = randn(NumberOfAtoms,1).* 20e-6;
    case 2
        ProbDF  = GaussianDistribution(0, (5e-6)^2, positions);
        initialPositions = drawSamplesFromDistribution(NumberOfAtoms, positions, ProbDF); 
    case 3 
        ProbDF  = HigherOrderGaussianDistribution(20e-6, (12e-6)^2, 2, positions);
        initialPositions = drawSamplesFromDistribution(NumberOfAtoms, positions, ProbDF); 
    case 4 
        ProbDF  = MirroredFlatTopGaussianDistribution(20e-6, (12e-6)^2, positions);
        initialPositions = drawSamplesFromDistribution(NumberOfAtoms, positions, ProbDF); 
end        

figure(1)
clf
NumberOfBins = 100;
histogram(initialPositions*1e6,NumberOfBins,'DisplayName','Sampled')
hold on
plot(positions*1e6, ProbDF*NumberOfAtoms*numel(ProbDF)/NumberOfBins,'LineWidth',3,'DisplayName','Predicted')
xlim([-round(max(initialPositions(:))*1e6,1), round(max(initialPositions(:))*1e6,1)])
sgtitle(PositionDistributionTypes{ChoiceOfType}, 'FontSize', 14)
%sgtitle('(Mirrored) Flat Top Gaussian', 'FontSize', 14)
legend('FontSize', 14)

max_vel = 25e-3;
velocities = 0:(max_vel/(NumberOfAtoms-1)):max_vel;

VelocityDistributionTypes = {'Zero', 'InBuiltNormal', 'Uniform', 'Gaussian', 'Maxwell-Boltzmann'};
ChoiceOfType = 5;

switch ChoiceOfType
    case 1
        initialVelocities  = 0;
    case 2
        initialVelocities  = randn(NumberOfAtoms,1).* 20e-6;
    case 3
        ProbDF     = UniformDistribution(0, 8e-03, velocities);
        initialVelocities = drawSamplesFromDistribution(NumberOfAtoms, velocities, ProbDF);
    case 4 
        ProbDF     = GaussianDistribution(10e-4, 2e-7, velocities);
        initialVelocities = drawSamplesFromDistribution(NumberOfAtoms, velocities, ProbDF);
    case 5 
        ProbDF     = MaxwellBoltzmannDistribution(initialTemperature, velocities);
        initialVelocities = drawSamplesFromDistribution(NumberOfAtoms, velocities, ProbDF);
end

signflips = (-1) .* (rand(length(initialVelocities),1) > 0.5);
signflips(signflips == 0) = 1;
initialVelocities = initialVelocities .* signflips;

figure(2)
clf
NumberOfBins = 100;
histogram(initialVelocities*1e3,NumberOfBins,'DisplayName','Sampled')
PDFForDisp = ProbDF*NumberOfAtoms*numel(ProbDF)/NumberOfBins;
PDFForDisp = horzcat(flip(PDFForDisp), PDFForDisp);
VelocityRange = horzcat(-flip(velocities*1e3), velocities*1e3);
hold on
plot(VelocityRange, PDFForDisp,'LineWidth',3,'DisplayName','Predicted')
%xlim([0, round(max(initialVelocities(:))*1e3,1)])
sgtitle(['Maxwell-Boltzmann for T = ' num2str(initialTemperature*1e9) ' nK'], 'FontSize', 14)
%sgtitle('(Mirrored) Flat Top Gaussian', 'FontSize', 14)
legend('FontSize', 14)

%%
if ~DebugMode
    tRes        = 1e-6;                                 % Resolution for the ODE solver (s)
    t0          = 0;                                    % Starting time (s)
    WaitTime    = 0.5e-3;                               % Final time (s)
    tNumPoints  = floor(WaitTime/tRes)+1;               % Number of sample points in time between t0 and tf
    tspan = linspace(t0,WaitTime,tNumPoints);           % Solver calculates atom position for each of these timesteps in this time array
    CompressionCycles = 2;                              % CompressionCycles
    
    AllPositionsOverTime = zeros(CompressionCycles * tNumPoints, NumberOfAtoms);
    AllVelocitiesOverTime = zeros(CompressionCycles * tNumPoints, NumberOfAtoms);
    for idx = 1:CompressionCycles
        [CurrentPositionDistribution,CurrentVelocityDistribution] = Solver(tRes, t0, WaitTime, Trap, InitialPositionDistribution, InitialVelocityDistribution);
        AllPositionsOverTime((idx-1)*size(CurrentPositionDistribution,1)+1:(idx*size(CurrentPositionDistribution,1)),:) = CurrentPositionDistribution;  
        AllVelocitiesOverTime((idx-1)*size(CurrentVelocityDistribution,1)+1:(idx*size(CurrentVelocityDistribution,1)),:) = CurrentVelocityDistribution;  
        InitialPositionDistribution = CurrentPositionDistribution(end,:);
        InitialVelocityDistribution = InitialVelocityDistribution(randperm(numel(InitialVelocityDistribution)));
    end
end
%% Plotting
if DebugMode
    PlotPotential                                  = true;
    PlotSampling                                   = true;
    PlotTrajectories                               = false;
    PlotEnvelopesAndOptimalWaitingTimes            = false;
    PlotPhaseSpaceEvolution                        = false;
    PlotRMSpread                                   = false;
    PlotPhaseSpaceAnalysis                         = false;
else
    PlotPotential                                  = false;
    PlotSampling                                   = false;
    PlotTrajectories                               = false;
    PlotEnvelopesAndOptimalWaitingTimes            = false;
    PlotPhaseSpaceEvolution                        = true;
    PlotRMSpread                                   = true;
    PlotPhaseSpaceAnalysis                         = false;
end

if PlotPotential
    plotPotential(Wavelength, Trap) 
end
if PlotSampling
    plotSampling(NumberOfAtoms, initialTemperature, velocities)
end
if PlotTrajectories
    TimeForFullCompression = plotTrajectories(NumberOfAtoms, Trap, tspan, AllPositionsOverTime);
    %plotQuarterPeriods(Trap, initialPositions, TimeForFullCompression)
end   
if PlotEnvelopesAndOptimalWaitingTimes
    plotEnvelopesAndOptimalWaitingTimes(tspan, Xres, 50)
end
if PlotPhaseSpaceEvolution
    NumberOfBins = 100;
    plotPhaseSpaceEvolution(NumberOfAtoms, Trap, GroundStatePopulation, FinalTrapDepthInUnitsOfRecoilEnergy, initialTemperature, t0, tNumPoints, WaitTime, CompressionCycles, AllPositionsOverTime, AllVelocitiesOverTime, NumberOfBins)
end
if PlotRMSpread
   plotRMSSpreadEvolution(t0, tNumPoints, WaitTime, AllPositionsOverTime, NumberOfAtoms, FractionOfInitialPotential, GroundStatePopulation, initialTemperature, CompressionCycles)
end 

if PlotPhaseSpaceAnalysis
    MaxHoldTime = 3; %in ms
    plotPhaseSpaceAnalysis(Trap, NumberOfAtoms, tspan, Xres, Vres, initialTemperature, MaxHoldTime)
end
%% Plotting functions
function plotPotential(Wavelength, Trap)
    PhysicsConstants;
    x = -200*Wavelength:1e-2*Wavelength:200*Wavelength;
    figure(1)
    clf
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    set(gcf, 'defaultAxesColorOrder', [[0 0 0];[0 0 0]]);
    Potential = @(x,Trap) Trap.U0 .* exp(-2 *((x/Trap.w0).^2));
    %plot(x*1e6,-cumsum(force(x,Trap)), 'Color', colours{2}, 'LineWidth', 2) % plot the potential
    plot(x*1e6, (Potential(x,Trap)./BoltzmannConstant) .*1e6, 'Color', colours{2}, 'LineWidth', 2) % plot the potential
    xlabel('Position (um)','FontSize', 14)
    ylabel('Potential (uK)','FontSize', 14)
    xlim([min(x) max(x)]*1e6)
    yyaxis right
    plot(x*1e6,force(x,Trap)*1e23, 'Color', colours{4}, 'LineWidth', 2) % plot the potential
    ylabel('Force (x 10^{-23} N)','FontSize', 14)
    legend({'Potential', 'Force'},'FontSize', 14)
    sgtitle(['Vertical Dipole Trap of depth = ' num2str(abs(Trap.U0InTemperature*1e6),'%.2f') ' uK'])
end
function plotSampling(NumberOfAtoms, initialTemperature, velocities)
% - plot Distribution
figure(2)

clf
ProbDFArray = cell(1,3);

ProbDFArray{1} = MaxwellBoltzmannDistribution(initialTemperature, velocities);
ProbDFArray{2} = GaussianDistribution(0, 2e-7, velocities);
ProbDFArray{3} = UniformDistribution(0, 8e-03, velocities);
ProbDFName = {'MaxwellBoltzmann','Gaussian','Uniform'};

for kk = 1:3
    ProbDF = ProbDFArray{kk};
    subplot(1,3,kk)
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    set(gcf, 'defaultAxesColorOrder', [[0 0 0];[0 0 0]]);
    initialVelocities = drawSamplesFromDistribution(NumberOfAtoms, velocities, ProbDF);
    signflips = (-1) .* (rand(length(initialVelocities),1) > 0.5);
    signflips(signflips == 0) = 1;
    initialVelocities = initialVelocities .* signflips;
    NumberOfBins = 100;
    histogram(initialVelocities*1e3,NumberOfBins, 'FaceAlpha', 0.1,'DisplayName','Sampled')
    xlim([round(min(initialVelocities(:))*1e3,1), round(max(initialVelocities(:))*1e3,1)])
    ylabel('Probability','FontSize', 14)
    title(ProbDFName{kk},'FontSize', 14)
    legend({'Sampled Velocities', 'PDF'})
end

sgtitle('Sampling from a Distribution','FontSize', 18)
end
function TimeForFullCompression = plotTrajectories(NumberOfAtoms, Trap, tspan, Xres)
% Plot all trajectories
figure(3)
clf
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
TimeForFullCompression = zeros(NumberOfAtoms,1);
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
for Index = 1:NumberOfAtoms
    plot(tspan*1e3, Xres(:,Index).*1e6, 'HandleVisibility', 'Off');
    hold on
    %ZC = findAllZeroCrossings(tspan*1e3,Xres(:,Index).*1e6);
    %TimeForFullCompression(Index) = ZC(1);
end
%MeanTime = mean(TimeForFullCompression);
%line([min(TimeForFullCompression) min(TimeForFullCompression)],[-10 10],'Color',colours{1},'LineStyle','--')
%line([MeanTime MeanTime],[-10 10],'Color',colours{1},'LineWidth',1.5)
%line([max(TimeForFullCompression) max(TimeForFullCompression)],[-10 10],'Color',colours{1},'LineStyle','--')
clear Index
sgtitle(['Trajectory of ' num2str(NumberOfAtoms) ' atoms at different starting positions in the VDT of depth: ' num2str(abs(Trap.U0InTemperature*1e6),'%.2f') ' uK']);
ylabel('Position (um)','FontSize', 14)
xlabel('Time (ms)','FontSize', 14)
%legend({['Minimum quarter period (' num2str(min(TimeForFullCompression),'%.3f') ' ms)'], ['Mean quarter period (' num2str(MeanTime,'%.3f') ' ms)'], ['Maximum quarter period (' num2str(max(TimeForFullCompression),'%.3f') ' ms)']}, 'FontSize', 14)
grid on
hold off
sgtitle(['Trajectory of ' num2str(NumberOfAtoms) ' atoms at different starting positions in the VDT of depth: ' num2str(abs(Trap.U0InTemperature*1e6),'%.2f') ' uK']);
end
function plotQuarterPeriods(Trap, initialPositions, TimeForFullCompression)
figure(4)
clf
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
set(gcf, 'defaultAxesColorOrder', [colours{1};colours{5}]);
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
%yyaxis left
[initialPositions, sortIdx] = sort(initialPositions, 'ascend');
TimeForFullCompression = TimeForFullCompression(sortIdx);
plot(initialPositions*1e6, TimeForFullCompression, '--o', 'Color',colours{1})
sgtitle(['For atoms starting at different positions in the VDT of depth: ' num2str(abs(Trap.U0InTemperature*1e6),'%.2f') ' uK']);
xlabel('Starting Position (um)','FontSize', 14)
ylabel('Quarter period [ms]','FontSize', 14)
grid on
% yyaxis right
% plot(initialPositions*1e6, 1E3./(4.*TimeForFullCompression), '--o', 'Color',colours{5})
% ylabel('Trapping Frequency [Hz]','FontSize', 14)
% legend({'Quarter period','Trapping Frequency'}, 'FontSize', 14)
end
function plotEnvelopesAndOptimalWaitingTimes(tspan, Xres, UptoInitPos)
    %%
    figure(5);
    clf;
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    initialPositions = Xres(1,:);
    [initialPositions, sortIdx] = sort(initialPositions, 'ascend');
    Xres = Xres(:,sortIdx);
    envelopes = zeros(size(Xres));
    optimalwaitingtimes  = zeros(size(Xres,2)-1,1);
    [~,OriginIndex]=min(abs(round(initialPositions.*1e6)));
    for idx = 1:size(Xres,2)-1
        FirstPosition = round(Xres(1,idx).*1e6);
        if  FirstPosition > 0
            Limit = find(round(initialPositions.*1e6)==FirstPosition, 1, 'first');
            envelopes(:,idx) = max(abs(Xres(:,OriginIndex:Limit)),[],2).*1e6;
        elseif FirstPosition < 0
            Limit = find(round(initialPositions.*1e6)==FirstPosition, 1, 'first');
            envelopes(:,idx) = max(abs(Xres(:,Limit:OriginIndex)),[],2).*1e6;
        end
        [~, index] = min(envelopes(:,idx));
        optimalwaitingtimes(idx) = tspan(index)*1e3;
    end
    idx = find(round(initialPositions.*1e6)==UptoInitPos, 1, 'first');
    idx = idx(1);
    plot(tspan*1e3, envelopes(:,idx), 'LineWidth', 5, 'Color',colours{2});
    hold on
    plot(tspan*1e3, abs(Xres).*1e6,'HandleVisibility', 'Off');
    xlabel('Time (ms)','FontSize', 14)
    ylabel('Position (um)','FontSize', 14)
    legend({['Upto ' num2str(initialPositions(idx)*1e6) ' um: \tau_{wait}^{opt} ~ ' num2str(optimalwaitingtimes(idx), '%.2f') ' ms']}, 'FontSize', 14)
    sgtitle('Optimal waiting times for different capture ranges');
    %%
    figure(6);
    clf
    plot(initialPositions(2:end-1)*1e6,optimalwaitingtimes(1:size(Xres, 2)-2), '--o', 'LineWidth', 2)
    xlabel('Capture Range (um)','FontSize', 14)
    ylabel('Optimal Waiting Time (ms)','FontSize', 14)
    sgtitle('Optimal waiting times for different capture ranges');
end
function plotPhaseSpaceEvolution(NumberOfAtoms, Trap, GroundStatePopulation, FinalTrapDepthInUnitsOfRecoilEnergy, initialTemperature, t0, tNumPoints, WaitTime, CompressionCycles, Xres, Vres, NumberOfBins)
    figure(7)
    clf
    PhysicsConstants;
    set(gcf, 'Units', 'normalized');
    set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    F = struct('cdata',[],'colormap',[]);
    frame = 0;
    Potential = @(x,Trap) Trap.U0 .* exp(-2 *((x/Trap.w0).^2));
    Pos = linspace(-max(Xres(:)), max(Xres(:)), size(Xres,2));
    EscapeVelocities = sqrt(2.*abs(Potential(Pos, Trap))./Cs133Mass)*1e3;
    tNumPoints   = CompressionCycles * tNumPoints;
    tspan        = linspace(t0,CompressionCycles * WaitTime,tNumPoints);
    for Time = 1:30:tNumPoints
        positions  = Xres(Time,:).*1e6;
        velocities = Vres(Time,:).*1e3;
        sb1 = subplot(4,4,[1,9]);
        h1 = histogram(velocities(:), NumberOfBins,'Normalization','probability');
        h1.FaceAlpha =  0.1;
        h1.FaceColor = colours{2};
        ylabel('Velocity Distribution','FontSize', 12)
        [MaxValue,~] = max(h1.Values);
%         hold on
%         [pdfey,pdfex] = ksdensity(velocities(:));
%         plot(pdfex, (pdfey./max(pdfey)).*MaxValue, 'Color', colours{6}, 'LineStyle', '--');
        set(sb1, 'Box', 'off', 'Color', 'none')
        set(sb1, 'xlim', [-round(max(Vres(:))*1e3, 1) round(max(Vres(:))*1e3, 1)])
        set(sb1, 'ylim', [0 MaxValue])
        set(sb1, 'XDir', 'reverse')
        camroll(sb1,90)
        hold off
        sb2 = subplot(4,4,[2.15,12]);
        plot2DHistogram(positions(:),  velocities(:), ...
            'nbins', NumberOfBins, ...
            'PositionLimits', [-round(max(Xres(:))*1e6, 1) round(max(Xres(:))*1e6, 1)],...
            'VelocityLimits', [-round(max(Vres(:))*1e3, 1) round(max(Vres(:))*1e3, 1)],...
            'CountDensity', false);
        colorbar
        hold on
        plot(Pos*1e6, EscapeVelocities, 'Color', [1 1 1], 'LineStyle', '--')
        plot(Pos*1e6, -EscapeVelocities, 'Color', [1 1 1], 'LineStyle', '--')
        text(round(min(Xres(:))*1e6, 1)+2, -round(max(Vres(:))*1e3, 1)+5, ['Time : ' num2str(round(tspan(Time)*1e3, 1)) ' ms'], 'Color', [1 1 1],'FontSize', 14)
        for idx = 1:CompressionCycles
            if tspan(Time) >= (idx-1)*WaitTime && tspan(Time) <= idx*WaitTime
                text(round(max(Xres(:))*1e6, 1)-8, -round(max(Vres(:))*1e3, 1)+5, ['Run : ' num2str(idx)], 'Color', [1 1 1],'FontSize', 14)
            end
        end
        xlabel('Position (\mum)','FontSize', 14)
        ylabel('Velocity (mm/s)','FontSize', 14)
        xlim([-round(max(Xres(:))*1e6, 1) round(max(Xres(:))*1e6, 1)])
        ylim([-round(max(Vres(:))*1e3, 1) round(max(Vres(:))*1e3, 1)])
        grid on
        sb3 = subplot(4,4,[14.15,15.7], 'color', 'none');
        h2 = histogram(positions(:), NumberOfBins,'Normalization','probability');
        h2.FaceAlpha =  0.1;
        h2.FaceColor = colours{2};
        ylabel('Position Distribution','FontSize', 12)
        [MaxValue,~] = max(h2.Values);
%         hold on
%         [pdfey,pdfex] = ksdensity(positions(:));
%         plot(pdfex, (pdfey./max(pdfey)).*MaxValue, 'Color', colours{6}, 'LineStyle', '--');
        set(sb3, 'Box', 'off', 'Color', 'none')
        set(sb3, 'xlim', [-round(max(Xres(:))*1e6, 1) round(max(Xres(:))*1e6, 1)])
        set(sb3, 'ylim', [0 MaxValue])
        hold off
        sgtitle(['Ground State Pop. = ' num2str(GroundStatePopulation) ' --> (Uniform) Initial Temperature = ' num2str(initialTemperature*1e6, '%.1f') ' \muK; Release Depth = ' num2str(FinalTrapDepthInUnitsOfRecoilEnergy) ' E_R']);
        frame = frame+1;
        F (frame) = getframe (gcf);
        drawnow
    end
    clear Index
    hold off
    writerObj = VideoWriter (['PhaseEvolution_N' num2str(NumberOfAtoms) '.avi']);
    writerObj.FrameRate = 10;
    writerObj.Quality = 100;
    open (writerObj);
    for i = 1: length (F)
        frame = F (i);
        writeVideo (writerObj, frame);
    end
    % close the writer object
    close (writerObj);
end
function plotRMSSpreadEvolution(t0, tNumPoints, WaitTime, Xres, NumberOfAtoms, FractionOfInitialPotential, GroundStatePopulation, initialTemperature, CompressionCycles)
    figure(8)
    clf
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    set(gcf, 'defaultAxesColorOrder', [colours{1};colours{5}]);
    set(gcf, 'Units', 'normalized');
    set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
    RMSSpread = zeros(size(Xres, 1),1);
    tNumPoints   = CompressionCycles * tNumPoints;
    tspan        = linspace(t0,CompressionCycles * WaitTime,tNumPoints);
    for idx = 1:size(Xres, 1)
        RMSSpread(idx) = rms(Xres(idx,:));
    end
    plot(tspan*1e3, RMSSpread.*1e6, 'Color',colours{6})
    ylim([0 max(RMSSpread)*1e6])
    text(1, min(RMSSpread)*1e6 -10, ['Number of Atoms:' num2str(NumberOfAtoms)],'FontSize', 14)
    text(1, min(RMSSpread)*1e6 -12, ['GS Population:' num2str(100*GroundStatePopulation) '%'],'FontSize', 14)
    text(1, min(RMSSpread)*1e6 -14, ['Initial Temperature (AFTER adiabatic release):' num2str(initialTemperature*1e9) ' nK'],'FontSize', 14)
    sgtitle(['Time evolution of RMS spread of atoms']);
    xlabel('Time (ms)','FontSize', 14)
    ylabel('RMS Spread [\mum]','FontSize', 14)
    legend({['Relative Release Depth:' num2str(100*FractionOfInitialPotential) '%']},'FontSize', 14)
    grid on
end
function plotRMSSpreadEvolutionForDifferentTemps(tspan, Spreads, NumberOfAtoms, FractionOfInitialPotential, GroundStatePopulation)
    figure(9)
    clf
    %colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    %set(gcf, 'defaultAxesColorOrder', [colours{1};colours{5}]);
    set(gcf, 'Units', 'normalized');
    set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
    for Index = 1:size(Spreads, 2)
        plot(tspan*1e3, Spreads{Index}.*1e6)
        hold on
    end
    ylim([0 max(max(cell2mat(Spreads)))*1e6])
    sgtitle(['Time evolution of RMS spread of atoms']);
    xlabel('Time (ms)','FontSize', 14)
    ylabel('RMS Spread [\mum]','FontSize', 14)
    legend(strcat(['#:' num2str(NumberOfAtoms) '; P_{GS}:'], strsplit(num2str(100.*GroundStatePopulation)), '%; ', ['U_{release}:' num2str(100*FractionOfInitialPotential) '%']),'FontSize', 14)
    grid on
end
function plotRMSSpreadEvolutionForDifferentReleaseDepths(tspan, Spreads, NumberOfAtoms, FractionOfInitialPotential, GroundStatePopulation)
    figure(9)
    clf
    %colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    %set(gcf, 'defaultAxesColorOrder', [colours{1};colours{5}]);
    set(gcf, 'Units', 'normalized');
    set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
    for Index = 1:size(Spreads, 2)
        plot(tspan*1e3, Spreads{Index}.*1e6)
        hold on
    end
    ylim([0 max(max(cell2mat(Spreads)))*1e6])
    sgtitle(['Time evolution of RMS spread of atoms']);
    xlabel('Time (ms)','FontSize', 14)
    ylabel('RMS Spread [\mum]','FontSize', 14)
    legend(strcat(['#:' num2str(NumberOfAtoms) '; P_{GS}:' num2str(100.*GroundStatePopulation) '%; U_{release}:'], strsplit(num2str(100.*FractionOfInitialPotential)), '%'),'FontSize', 14)
    grid on
end
function plotConvolvedDistributionSampling(NumberOfAtoms, initialTemperature, velocities)
% - plot Distribution
figure(9)
clf
PhysicsConstants;
InitialLongitudinalTrapFrequency = 56; %in kHz
FractionOfInitialPotential = 0.05;
InitialTrapDepth = LatticeProperties.estimateTrapDepthFromHeatingSidebandFreqDT1DT3(InitialLongitudinalTrapFrequency);
TrapDepth = sqrt(FractionOfInitialPotential) * InitialTrapDepth;
results = LatticeProperties.estimateTrapFreq33FromTrapDepth(0,0,0,0,TrapDepth);
FinalLongitudinalTrapFrequency = results(1);
DeltaPQuantumMechanical = sqrt(0.5*Cs133Mass*PlanckConstantReduced*FinalLongitudinalTrapFrequency)./Cs133Mass;
ProbDFArray = cell(1,3);
ProbDFArray{1}    = MaxwellBoltzmannDistribution(initialTemperature, velocities);
ProbDFArray{2}    = GaussianDistribution(0, DeltaPQuantumMechanical, velocities);
ProbDFArray{3}    = conv(ProbDFArray{1}, ProbDFArray{2});
ProbDFName = {'MaxwellBoltzmann','Gaussian','Convolution'};
for kk = 1:3
    ProbDF = ProbDFArray{kk};
    MaxVelocity = max(velocities);
    subplot(1,3,kk)
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    set(gcf, 'defaultAxesColorOrder', [[0 0 0];[0 0 0]]);
    ProbDF     = ProbDF./sum(ProbDF);
    if kk==3
        initialVelocities = drawSamplesFromDistribution(NumberOfAtoms, 0:MaxVelocity/(length(ProbDF)-1):MaxVelocity, ProbDF);
        histogram(initialVelocities.*1e6, 'FaceAlpha', 0.1)
        xlabel('Velocities(m/s)','FontSize', 14)
        ylabel('Frequency','FontSize', 14)
        legend({'Sampled Velocities', 'PDF'})
        yyaxis right
    end
    stairs((0:(MaxVelocity/(length(ProbDF)-1)):MaxVelocity).*1e6, ProbDF, 'Color', colours{4}, 'LineWidth', 2)
    xlabel('Velocities(um/s)','FontSize', 14)
    ylabel('Probability','FontSize', 14)
    if kk~=3
        yyaxis right
        plot((0:(MaxVelocity/(length(ProbDF)-1)):MaxVelocity).*1e6, cumsum(ProbDF))
    end
    title(ProbDFName{kk},'FontSize', 14)
end
sgtitle('Sampling from a Distribution','FontSize', 18)
end
function plotPhaseSpaceAnalysis(Trap, NumberOfAtoms, tspan, Xres, Vres, initialTemperature, MaxHoldTime)
    PhysicsConstants;
    figure(10)
    clf
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    set(gcf, 'defaultAxesColorOrder', [colours{1};colours{5}]);
    set(gcf, 'Units', 'normalized');
    set(gcf, 'OuterPosition', [0.025 0.1 0.9 0.9]);
    
    PositionSpreads = zeros(size(Xres, 1),1);
    for idx = 1:size(Xres, 1)
        PositionSpreads(idx) = rms(Xres(idx,:));
    end
    AtomCounts = zeros(size(Xres, 1),1);
    for Idx = 1:length(PositionSpreads)
       AtomCounts(Idx) = sum(abs(Xres(Idx, :)) <= PositionSpreads(Idx));
    end    
    PositionMeans = zeros(size(Xres, 1),1);
    for idx = 1:size(Xres, 1)
        PositionMeans(idx) = mean(Xres(idx,:));
    end
    VelocitySpreads = zeros(size(Vres, 1),1);
    for idx = 1:size(Xres, 1)
        VelocitySpreads(idx) = rms(Vres(idx,:));
    end
    VelocityMeans = zeros(size(Vres, 1),1);
    for idx = 1:size(Xres, 1)
        VelocityMeans(idx) = mean(Vres(idx,:));
    end
    Periods = tspan.*1e3;
    [~,Idx]=min(abs(Periods-MaxHoldTime));
    Periods = Periods(1:Idx);
    PositionSpreads = PositionSpreads(1:Idx);
    AtomCounts = AtomCounts(1:Idx);
    VelocitySpreads = VelocitySpreads(1:Idx);
    PositionMeans = PositionMeans(1:Idx);
    VelocityMeans = VelocityMeans(1:Idx);
    
    %'Amplitude','Frequency','Phase','Offset'
    
    %Sin
    Sin_Func_Handle = @(p,x)  (p(4)+ (p(1) * sin((2*pi*p(2)*x)+p(3)))); 
    
    %SinSquared
    SinSq_Func_Handle = @(p,x)  (p(4)+ (p(1) * 1/2 * (1-cos((2*pi*p(2)*x)+p(3))))); 
    
    subplot(2,3,1)
    plot(Periods, PositionSpreads.*1e6, 'Color',colours{1}, 'LineWidth', 2)
    [~,Idx]=min(abs(PositionSpreads));
    line([Periods(Idx) Periods(Idx)], [0 max(PositionSpreads.*1e6)], 'Color',colours{1}, 'LineStyle', '--', 'LineWidth', 2);
    text(Periods(Idx)-0.1, 0.8*(max(PositionSpreads)-min(PositionSpreads))*1e6, ['t = ' num2str(Periods(Idx)) ' ms'], 'Rotation', 90, 'FontSize', 14)
    ylim([0 max(PositionSpreads)*1e6+1])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('\mum','FontSize', 14)
    legend({'x_{RMS}'},'FontSize', 14)
    grid on
    title('RMS of positions over time')
    
    subplot(2,3,2)
    plot(Periods, VelocitySpreads.*1e3, 'Color',colours{2}, 'LineWidth', 2)
    VRMSClassical = 1e3 * sqrt((BoltzmannConstant*initialTemperature/Cs133Mass));
    line([0 max(Periods)], [VRMSClassical VRMSClassical], 'LineStyle', '--', 'LineWidth', 2)
    text(max(Periods)-2.5, VRMSClassical+9, ['v^{t=0}_{RMS} = ' num2str(VelocitySpreads(1)*1e3) ' mm/s'], 'FontSize', 14)
    text(max(Periods)-2.5, VRMSClassical+3, ['Expected v^{t=0}_{RMS} = ' num2str(VRMSClassical) ' mm/s'], 'FontSize', 14)
    ylim([0 max(VelocitySpreads.*1e3)+5])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('mm/s','FontSize', 14)
    legend({'v_{RMS}'},'FontSize', 14)
    grid on
    title('RMS of velocities over time')
    
%     subplot(2,3,3)
%     grid on
%     title('Cross-Check')
    
    subplot(2,3,3)
    plot(Periods, PositionMeans.*1e6, 'Color',colours{4},  'LineWidth', 2)
    hold on
    %plot(Periods, Sin_Func_Handle([max(PositionMeans).*1e6, Trap.LongitudinalTrapFrequencyinHz, pi/2, 0], Periods), 'Color',colours{1}, 'LineWidth', 2, 'LineStyle', '--');
    ZC  = findAllZeroCrossings(Periods,PositionMeans);
    for Idx = 1:length(ZC)
        line([ZC(Idx) ZC(Idx)], [min(PositionMeans.*1e6) max(PositionMeans.*1e6)], 'Color',colours{4}, 'LineStyle', '--', 'LineWidth', 2);
        text(ZC(Idx)-0.1, min(PositionMeans)*1e6, ['t = ' num2str(ZC(Idx)) ' ms'], 'Rotation', 90, 'FontSize', 14)
    end
    xlim([0 MaxHoldTime])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('\mum','FontSize', 14)
    %legend({'x_{Mean}', 'Expected'},'FontSize', 14)
    grid on
    title('Mean of positions over time')
    
    subplot(2,3,4)
    plot(Periods, VelocityMeans.*1e3, 'Color',colours{5}, 'LineWidth', 2)
    hold on
    %plot(Periods, Sin_Func_Handle([(max(PositionMeans) * Trap.LongitudinalTrapFrequencyinHz), Trap.LongitudinalTrapFrequencyinHz, pi, 0], Periods), 'Color',colours{1}, 'LineWidth', 2, 'LineStyle', '--');
    for Idx = 1:length(ZC)
        line([ZC(Idx) ZC(Idx)], [min(VelocityMeans.*1e3) max(VelocityMeans.*1e3)], 'Color',colours{5}, 'LineStyle', '--', 'LineWidth', 2);
    end
    xlim([0 MaxHoldTime])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('mm/s','FontSize', 14)
    %legend({'v_{Mean}', 'Expected'},'FontSize', 14)
    grid on
    title('Mean of velocities over time')
    
    subplot(2,3,5)
    NormalizedAtomCounts = AtomCounts./NumberOfAtoms; 
    plot(Periods, NormalizedAtomCounts, 'Color',colours{6}, 'LineWidth', 2)
    for Idx = 1:length(ZC)
        line([ZC(Idx) ZC(Idx)], [min(NormalizedAtomCounts) max(NormalizedAtomCounts)], 'Color',colours{6}, 'LineStyle', '--', 'LineWidth', 2);
    end
    xlim([0 MaxHoldTime])
    ylim([0 1])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('Normalized Atom Count','FontSize', 14)
    grid on
    title('# of Atoms within the position RMS')
    
    sgtitle(['Quantitative analysis of phase space dynamics']);
end
function plot2DHistogram(x,y,varargin)
% hist2d(x,y,varargin) this plots the phase-space-density for ... 
    
    p = inputParser;
    addRequired(p,  'Positions',            @isnumeric)
    addRequired(p,  'Velocities',           @isnumeric)
    addParameter(p, 'nbins',           10,  @isscalar)
    addParameter(p, 'PositionLimits',  [-60 60],  @(x) isDataArray(x) && isnumeric(x))
    addParameter(p, 'VelocityLimits',  [-80 80],  @(x) isDataArray(x) && isnumeric(x))
    addParameter(p, 'CountDensity', false,  @islogical)
    parse(p,x,y,varargin{:})
    
    x = p.Results.Positions;
    y = p.Results.Velocities;
    nbins = p.Results.nbins;
    plim  = p.Results.PositionLimits;
    vlim  = p.Results.VelocityLimits;
    cdflag  = p.Results.CountDensity;
    
    xedges = linspace(min(plim),max(plim),nbins+1);
    yedges = linspace(min(vlim),max(vlim),nbins+1);
    
    % computes bins vectors (as middle of each edges couples)
    xbins = mean(cat(1,xedges(1:end-1),xedges(2:end)));
    ybins = mean(cat(1,yedges(1:end-1),yedges(2:end)));

    % computes bins width vectors and area matrix
    xbw = diff(xedges);
    ybw = diff(yedges);
    [xx,yy] = meshgrid(xbw,ybw);
    a = xx.*yy;
    
    % initiate the result matrix
    n = zeros(length(ybins),length(xbins));
    % main loop to fill the matrix with element counts
    for i = 1:size(n,1)
        k = find(y >= yedges(i) & y < yedges(i+1));
        for j = 1:size(n,2)
            n(i,j) = length(find(x(k) >= xedges(j) & x(k) < xedges(j+1)));
        end
    end
    
    % normalize options
    if cdflag
        n = n./a; % count density
    end
    imagesc(xbins,ybins,n)
    %hold on
    %plot(x,y,'.k','MarkerSize',10)
    %hold off	
end
%% helper functions
function ret = odefcn(~, x, Trap)
% ret = odefcn(~, x, Trap) implements the ODE for the simulation
    PhysicsConstants;
    ret = zeros(2,1);
    ret(1) = x(2); % Velocity
    ret(2) = force(x(1), Trap)/Cs133Mass;
end % - ODE to solve
function ret = force(x, Trap, varargin)
% ret = force(x, Trap) calculates the force [N] that the dipole-potential 
% exerts on a CS133-atom at distance x [m] from its center 
    PhysicsConstants;
    narginchk(2,3)
    if nargin == 2
        Type = 'Gaussian';
    else
        Type = varargin{1};
    end
    switch Type
        case 'Parabolic'
            ret = - Cs133Mass * Trap.TransverseTrappingFrequency^2 * x;
        case 'Gaussian'
            ret = Trap.U0 * 4 * x/(Trap.w0.^2) .* exp(-2 *((x/Trap.w0).^2));
    end  
end     % - Dipole force given a Gaussian potential  
function ret = potential(x, Trap, varargin)
    PhysicsConstants;
    PhysicsConstants;
    narginchk(2,3)
    if nargin == 2
        Type = 'Gaussian';
    else
        Type = varargin{1};
    end
    switch Type
        case 'Parabolic'
            ret = 0.5 * Cs133Mass * Trap.TransverseTrappingFrequency^2 .* x.^2;
        case 'Gaussian'
            ret = Trap.U0 .* exp(-2 *((x/Trap.w0).^2));
    end   
end
function [Xres, Vres] = Solver(tRes, t0, tf, Trap, AtomPositionDistribution, AtomVelocityDistribution)
    tNumPoints  = floor(tf/tRes)+1;     % Number of sample points in time between t0 and tf
    tspan = linspace(t0,tf,tNumPoints); % Solver calculates atom position for each of these timesteps in this time array
    NumberOfAtoms = length(AtomPositionDistribution);
    Xres = zeros(length(tspan),NumberOfAtoms);
    Vres = zeros(length(tspan),NumberOfAtoms);  
    progressbar  = parforNotifications();
    progressbar.PB_start(NumberOfAtoms,'Message',['Computing evolution for ' num2str(NumberOfAtoms,'%.0f') ' atoms:']);
    parfor Index = 1:NumberOfAtoms
        InitialConditions = [AtomPositionDistribution(Index) AtomVelocityDistribution(Index)]; % [initial position, initial velocity]
        [res] = ode5(@(t, x) odefcn(t, x, Trap), tspan, InitialConditions);
        Xres(:,Index) = res(:,1);
        Vres(:,Index) = res(:,2);
        progressbar.PB_iterate();
    end
    clear Index
end
function ret = UniformDistribution(a, b, x)
    ret = zeros(length(x),1);
    ret(x>a & x<b) = 1/(b - a);
end
function ret = GaussianDistribution(mean, var, x)
    % return the gauss distribution with parameters mean, var.
    % p = gaussdis(mean, var, x);
    ret = 1/sqrt(2*pi*var)*exp(-((x-mean).^2)/(2*var));
    ret = ret./sum(ret);
end
function ret = HigherOrderGaussianDistribution(mean, var, exponent, x)
    ret = exp(-((x-mean).^2./(2*var)).^exponent);
    ret = ret./sum(ret);
end
function ret = MirroredFlatTopGaussianDistribution(mean, var, x)
    IndexOfPositionClosestToZero = find(x(:).*circshift(x(:), [-1 0]) <= 0, 1, 'first');
    PositionClosestToZero = x(IndexOfPositionClosestToZero);     
    if PositionClosestToZero < 0
        IndexOfPositionClosestToZero = IndexOfPositionClosestToZero + 1;
    end
    x = x(IndexOfPositionClosestToZero:end);
    ret = exp(-max((x-mean),0).^2./(2*var));
    ret = horzcat(flip(ret),ret);
    ret = ret./sum(ret);
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
function ZC  = findAllZeroCrossings(x,y)
% ZC = findAllZeroCrossings(x,y) finds all Zero-crossing of the function y = f(x)
%
% 
% Remark:
%   findAllZeroCrossings has been modified so as to avoid the following error "Error using griddedInterpolant
%   The grid vectors are not strictly monotonic increasing."
%
%   The following modifications are based off of this MatLab forum answer:
%   https://www.mathworks.com/matlabcentral/answers/283718-grid-vectors-not-strictly-monotonic-increasing
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns Approximate Zero-Crossing Indices Of Argument Vector
zxidx = zci(y);
if ~isempty(zxidx)
    for k1 = 1:numel(zxidx)
        idxrng = max([1 zxidx(k1)-1]):min([zxidx(k1)+1 numel(y)]);
        xrng = x(idxrng);
        yrng = y(idxrng);
        % Beginning of findAllZeroCrossings2 modifications. The naming conventions follow
        % those in the referenced MatLab forum, except that "X" is "yrng" and
        % "Y" is "xrng".
        [yrng2, ~, jyrng] = unique(yrng); %yrng is a new array containing the unique values of yrng. jyrng contains the indices in yrng that correspond to the original vector. yrng = yrng2(jyrng)
        xrng2 = accumarray(jyrng, xrng, [], @mean); %This function creates a new array "xrng2" by applying the function "@mean" to all elements in "xrng" that have identical indices in "jyrng". Any elements with identical X values will have identical indices in jyrng. Thus, this function creates a new array by averaging values with identical X values in the original array.
        ZC(k1) = interp1( yrng2(:), xrng2(:), 0, 'linear', 'extrap' );
    end
else
    warning('No zero crossings found!')
    ZC = nan;
end
end