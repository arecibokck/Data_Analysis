function plotSpreadEvolution(Solver)
    f_h = DQSIMhelper.getFigureByTag('SpreadTimeEvolution');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'Time evolution of RMS Spread';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 955 761];
    
    n = Solver.numberOfAtoms;
    GroundStatePopulation = Solver.InitialDistributionParameters.GroundStatePopulation;
    FractionOfInitialPotential = Solver.InitialDistributionParameters.FractionOfInitialPotential;
    initialTemperature = Solver.InitialDistributionParameters.initialTemperature;
    
    tspan = Solver.timeSpan;
    
    positions = Solver.simulationResults(:,:,1);
    RMSSpread = zeros(size(positions, 1),1);
    for idx = 1:size(positions, 1)
        RMSSpread(idx) = rms(positions(idx,:));
    end
    plot(tspan*1e3, RMSSpread.*1e6, 'Color', [0.6350, 0.0780, 0.1840])
    ylim([0 max(RMSSpread)*1e6+1])
    
    a_h = get(f_h, 'CurrentAxes');
    maxxlim = max(a_h.XLim);
    maxylim = max(a_h.YLim);
    text((maxxlim - (0.3 * maxxlim)) , (0.3 * maxylim), ['Number of Atoms:' num2str(n)],'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text((maxxlim - (0.28 * maxxlim)), (0.2 * maxylim), ['GS Population:' num2str(100*GroundStatePopulation) '%'],'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text((maxxlim - (0.7 * maxxlim)) , (0.1 * maxylim), ['Initial Temperature (AFTER adiabatic release):' num2str(initialTemperature*1e9) ' nK'],'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    
    sgtitle(['Time evolution of RMS spread of atoms']);
    xlabel('Time (ms)','FontSize', 14)
    ylabel('RMS Spread [\mum]','FontSize', 14)
    legend({['Relative Release Depth:' num2str(100*FractionOfInitialPotential) '%']},'FontSize', 14)
    grid on
end