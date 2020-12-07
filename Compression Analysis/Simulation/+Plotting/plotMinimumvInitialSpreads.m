function plotMinimumvInitialSpreads(Solver, InitialSpreads, MinimumSpreads, GroundStatePopulations)
    f_h = DQSIMhelper.getFigureByTag('MinimumSpreadsvInitialSpreads');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'MinimumSpreadsvInitialSpreads';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 760 600];
    
    n = Solver.numberOfAtoms;
    FractionOfInitialPotential= Solver.InitialDistributionParameters.FractionOfInitialPotential;
    
    Markers = {'o', '+', 'v', 's', 'p', 'h', 'd', '^'};
    for Index = 1:size(MinimumSpreads, 2)
        plot(InitialSpreads, MinimumSpreads{Index}.*1e6, strcat('--', Markers{Index}), 'MarkerSize', 5)
        hold on
    end
    xlabel('Initial Spread [\mum]','FontSize', 14)
    ylabel('First minimum of RMS Spread [\mum]','FontSize', 14)
    legend(strcat(['#:' num2str(n) '; P_{GS}:'], strsplit(num2str(100.*GroundStatePopulations)), '%; ', ['U_{release}:' num2str(100*FractionOfInitialPotential) '%']),'FontSize', 14)
    grid on
    sgtitle('First minimum of spread for different initial distributions');
end