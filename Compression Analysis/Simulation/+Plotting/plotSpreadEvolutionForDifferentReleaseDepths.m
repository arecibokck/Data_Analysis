function plotSpreadEvolutionForDifferentReleaseDepths(Solver, FractionsOfInitialPotential, Spreads)
    f_h = DQSIMhelper.getFigureByTag('SpreadsForDifferentReleaseDepths');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'Time evolution of RMS Spread for different release depths';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 760 600];
    
    GroundStatePopulation = Solver.InitialDistributionParameters.GroundStatePopulation;
    
    tspan = Solver.timeSpan;
    
    for Index = 1:size(Spreads, 2)
        plot(tspan*1e3, Spreads{Index}.*1e6)
        hold on
    end
    
    ylim([0 max(max(cell2mat(Spreads)))*1e6])
    sgtitle(['Time evolution of RMS spread of atoms']);
    xlabel('Time (ms)','FontSize', 14)
    ylabel('RMS Spread [\mum]','FontSize', 14)
    legend(strcat(['#:' num2str(NumberOfAtoms) '; P_{GS}:' num2str(100.*GroundStatePopulation) '%; U_{release}:'], strsplit(num2str(100.*FractionsOfInitialPotential)), '%'),'FontSize', 14)
    grid on
end
