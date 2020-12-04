function TimeForFullCompression = plotTrajectories(Solver)
    % Plot all trajectories
    f_h = DQSIMhelper.getFigureByTag('trajectories');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'Trajectories in Position Space';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 1.2493e+03 896];
    
    positions = Solver.simulationResults(:,:,1);
    n = Solver.numberOfAtoms;
    tspan = Solver.timeSpan;
    
    TimeForFullCompression = zeros(n,1);
    for Index = 1:n
        plot(tspan*1e3, positions(:,Index).*1e6, 'HandleVisibility', 'Off');
        hold on
        ZC = CompressionHelper.findAllZeroCrossings(tspan*1e3,positions(:,Index).*1e6);
        TimeForFullCompression(Index) = ZC(1);
    end
    
    MeanTime = mean(TimeForFullCompression);
    line([min(TimeForFullCompression) min(TimeForFullCompression)],[-10 10],'Color',[0, 0.4470, 0.7410],'LineStyle','--');
    line([MeanTime MeanTime],[-10 10],'Color',[0, 0.4470, 0.7410],'LineWidth',1.5);
    line([max(TimeForFullCompression) max(TimeForFullCompression)],[-10 10],'Color',[0, 0.4470, 0.7410],'LineStyle','--');
    clear Index
    sgtitle(['Trajectory of ' num2str(n) ' atoms at different starting positions in a dipole trap of depth: ' num2str(abs(Solver.potentialDepthInTemp),'%.2f') ' uK']);
    ylabel('Position (um)','FontSize', 14)
    xlabel('Time (ms)','FontSize', 14)
    legend({['Minimum quarter period (' num2str(min(TimeForFullCompression),'%.3f') ' ms)'], ['Mean quarter period (' num2str(MeanTime,'%.3f') ' ms)'], ['Maximum quarter period (' num2str(max(TimeForFullCompression),'%.3f') ' ms)']}, 'FontSize', 14)
    grid on
    hold off
    sgtitle(['Trajectory of ' num2str(n) ' atoms at different starting positions in a dipole trap of depth: ' num2str(abs(Solver.potentialDepthInTemp),'%.2f') ' uK']);
end