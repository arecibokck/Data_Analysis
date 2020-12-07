function plotQuarterPeriods(Solver, TimeForFullCompression)
    f_h = DQSIMhelper.getFigureByTag('QuarterPeriods');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children'))
        clf(f_h);
    end
    f_h.Name = 'Quarter periods for different starting positions';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 760 600];
    
    initialPositions  = Solver.initialPositions;
    
    [initialPositions, sortIdx] = sort(initialPositions, 'ascend');
    TimeForFullCompression = TimeForFullCompression(sortIdx);
    plot(initialPositions*1e6, TimeForFullCompression, '--o', 'Color',[0, 0.4470, 0.7410])
    sgtitle(['For atoms starting at different positions in the VDT of depth: ' num2str(abs(Solver.potentialDepthInTemp),'%.2f') ' uK']);
    xlabel('Starting Position (um)','FontSize', 14)
    ylabel('Quarter period [ms]','FontSize', 14)
    grid on
end
