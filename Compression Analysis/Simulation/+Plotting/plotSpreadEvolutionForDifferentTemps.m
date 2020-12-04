function plotSpreadEvolutionForDifferentTemps(Solver, GroundStatePopulations, Spreads)
    f_h = DQSIMhelper.getFigureByTag('SpreadsForDifferentTemps');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'Time evolution of RMS Spread for different temperatures';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 992 783];
    
    n = Solver.numberOfAtoms;
    FractionOfInitialPotential = Solver.InitialDistributionParameters.FractionOfInitialPotential;
    
    tspan = Solver.timeSpan;
    
    for Index = 1:size(Spreads, 2)
        plot(tspan*1e3, Spreads{Index}.*1e6)
        hold on
    end
    
    ylim([0 max(max(cell2mat(Spreads)))*1e6])
    sgtitle('Time evolution of RMS spread of atoms');
    xlabel('Time (ms)','FontSize', 14)
    ylabel('RMS Spread [\mum]','FontSize', 14)
    legend(strcat(['#:' num2str(n) '; P_{GS}:'], strsplit(num2str(100.*GroundStatePopulations)), '%; ', ['U_{release}:' num2str(100*FractionOfInitialPotential) '%']),'FontSize', 14, 'Location', 'southeast')
    grid on
end