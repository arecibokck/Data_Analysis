function plotMinimumSpreadPeriodvInitialPosition(Solver, initialPositions, MinimumSpreadPeriods)
    f_h = DQSIMhelper.getFigureByTag('MinimumSpreadsvInitialPositions');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'MinimumSpreadsvInitialPositions';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 760 600];
    
    GroundStatePopulation = Solver.InitialDistributionParameters.GroundStatePopulation;
    
    Markers = {'o', '+', 'v', 's', 'p', 'h', 'd', '^'};
    HarmonicQuarterPeriod = 0.25 * 1/(Trap.TrappingFrequency/(2*pi));
    for Index = 1:size(MinimumSpreadPeriods, 2)
        plot((initialPositions ./ Trap.w0), (MinimumSpreadPeriods{Index} ./ HarmonicQuarterPeriod), strcat('--', Markers{Index}), 'MarkerSize', 5)
        hold on
    end
    ylim([0 2])
    sgtitle('Variation of quarter periods for starting positions up to the beam waist');
    xlabel('x/w_o','FontSize', 14)
    ylabel('\tau / \tau_{Harmonic}','FontSize', 14)
    legend(strcat('P_{GS}:', strsplit(num2str(100.*GroundStatePopulation)), '% '),'FontSize', 14)
    grid on
end