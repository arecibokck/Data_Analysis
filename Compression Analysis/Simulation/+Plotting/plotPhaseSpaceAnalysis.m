function plotPhaseSpaceAnalysis(Solver,MaxHoldTime)
    f_h = DQSIMhelper.getFigureByTag('PhaseSpaceAnalysis');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'Analysis of Phase-Space Evolution';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 1584 865];
    
    n = Solver.numberOfAtoms;
    initialTemperature = Solver.InitialDistributionParameters.initialTemperature;
    
    tspan = Solver.timeSpan;
    pc = Solver.physicalConstants;
    
    positions  = Solver.simulationResults(:,:,1);
    velocities = Solver.simulationResults(:,:,2);
    
    PositionSpreads = zeros(size(positions, 1),1);
    for idx = 1:size(positions, 1)
        PositionSpreads(idx) = rms(positions(idx,:));
    end
    AtomCounts = zeros(size(positions, 1),1);
    for Idx = 1:length(PositionSpreads)
       AtomCounts(Idx) = sum(abs(positions(Idx, :)) <= PositionSpreads(Idx));
    end    
    PositionMeans = zeros(size(positions, 1),1);
    for idx = 1:size(positions, 1)
        PositionMeans(idx) = mean(positions(idx,:));
    end
    VelocitySpreads = zeros(size(velocities, 1),1);
    for idx = 1:size(positions, 1)
        VelocitySpreads(idx) = rms(velocities(idx,:));
    end
    VelocityMeans = zeros(size(velocities, 1),1);
    for idx = 1:size(positions, 1)
        VelocityMeans(idx) = mean(velocities(idx,:));
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
    plot(Periods, PositionSpreads.*1e6, 'Color',[0, 0.4470, 0.7410], 'LineWidth', 2)
    [~,Idx]=min(abs(PositionSpreads));
    line([Periods(Idx) Periods(Idx)], [0 max(PositionSpreads.*1e6)], 'Color',[0.8500, 0.3250, 0.0980], 'LineStyle', '--', 'LineWidth', 2);
    text(Periods(Idx)-0.1, 0.8*(max(PositionSpreads)-min(PositionSpreads))*1e6, ['t = ' num2str(Periods(Idx)) ' ms'], 'Rotation', 90, 'FontSize', 14)
    ylim([0 max(PositionSpreads)*1e6+1])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('\mum','FontSize', 14)
    legend({'x_{RMS}'},'FontSize', 14)
    grid on
    title('RMS of positions over time')
    
    subplot(2,3,2)
    plot(Periods, VelocitySpreads.*1e3, 'Color',[0.9290, 0.6940, 0.1250], 'LineWidth', 2)
    VRMSClassical = 1e3 * sqrt((pc.BoltzmannConstant*initialTemperature/pc.Cs133Mass));
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
    plot(Periods, PositionMeans.*1e6, 'Color',[0.4940, 0.1840, 0.5560],  'LineWidth', 2)
    hold on
    %plot(Periods, Sin_Func_Handle([max(PositionMeans).*1e6, Trap.LongitudinalTrapFrequencyinHz, pi/2, 0], Periods), 'Color',colours{1}, 'LineWidth', 2, 'LineStyle', '--');
    ZC  = CompressionHelper.findAllZeroCrossings(Periods,PositionMeans);
    for Idx = 1:length(ZC)
        line([ZC(Idx) ZC(Idx)], [min(PositionMeans.*1e6) max(PositionMeans.*1e6)], 'Color',[0.4940, 0.1840, 0.5560], 'LineStyle', '--', 'LineWidth', 2);
        text(ZC(Idx)-0.1, min(PositionMeans)*1e6, ['t = ' num2str(ZC(Idx)) ' ms'], 'Rotation', 90, 'FontSize', 14)
    end
    xlim([0 MaxHoldTime])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('\mum','FontSize', 14)
    %legend({'x_{Mean}', 'Expected'},'FontSize', 14)
    grid on
    title('Mean of positions over time')
    
    subplot(2,3,4)
    plot(Periods, VelocityMeans.*1e3, 'Color',[0.4660, 0.6740, 0.1880], 'LineWidth', 2)
    hold on
    %plot(Periods, Sin_Func_Handle([(max(PositionMeans) * Trap.LongitudinalTrapFrequencyinHz), Trap.LongitudinalTrapFrequencyinHz, pi, 0], Periods), 'Color',colours{1}, 'LineWidth', 2, 'LineStyle', '--');
    for Idx = 1:length(ZC)
        line([ZC(Idx) ZC(Idx)], [min(VelocityMeans.*1e3) max(VelocityMeans.*1e3)], 'Color',[0.4660, 0.6740, 0.1880], 'LineStyle', '--', 'LineWidth', 2);
    end
    xlim([0 MaxHoldTime])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('mm/s','FontSize', 14)
    %legend({'v_{Mean}', 'Expected'},'FontSize', 14)
    grid on
    title('Mean of velocities over time')
    
    subplot(2,3,5)
    NormalizedAtomCounts = AtomCounts./n; 
    plot(Periods, NormalizedAtomCounts, 'Color',[0.6350, 0.0780, 0.1840], 'LineWidth', 2)
    for Idx = 1:length(ZC)
        line([ZC(Idx) ZC(Idx)], [min(NormalizedAtomCounts) max(NormalizedAtomCounts)], 'Color',[0.6350, 0.0780, 0.1840], 'LineStyle', '--', 'LineWidth', 2);
    end
    xlim([0 MaxHoldTime])
    ylim([0 1])
    xlabel('t_{hold} (in ms)','FontSize', 14)
    ylabel('Normalized Atom Count','FontSize', 14)
    grid on
    title('# of Atoms within the position RMS')
    
    sgtitle(['Quantitative analysis of phase space dynamics']);
end

