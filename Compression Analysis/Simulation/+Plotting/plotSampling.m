function plotSampling(Solver, PositionDistribution, VelocityDistribution, NumberOfBins)
    % - plot Distribution
    f_h = DQSIMhelper.getFigureByTag('SamplingCheck');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'Sampling from a distribution';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 1.2947e+03 600];
    
    n = Solver.numberOfAtoms;
     
    initialPositions  = Solver.initialPositions;
    PositionPDFForDisp = PositionDistribution * n * numel(PositionDistribution)/NumberOfBins;
    
    initialVelocities = Solver.initialVelocities;
    VelocityPDFForDisp = VelocityDistribution * n * numel(VelocityDistribution)/NumberOfBins;
    
    TypeOfPositionDistribution = Solver.InitialDistributionParameters.TypeOfPositionDistribution;
    TypeOfVelocityDistribution = Solver.InitialDistributionParameters.TypeOfVelocityDistribution;
    
    initialTemperature = Solver.InitialDistributionParameters.initialTemperature;
    
    subplot(1,2,1)
    histogram(initialPositions*1e6,NumberOfBins, 'LineStyle', 'none', 'DisplayName','Sampled')
    hold on
    plot(Solver.positionGrid*1e6, PositionPDFForDisp, 'Color', [0.6350, 0.0780, 0.1840], 'LineStyle', '--', 'LineWidth', 1.5,'DisplayName','Predicted')
    xlim([-round(max(initialPositions(:))*1e6,1), round(max(initialPositions(:))*1e6,1)])
    xlabel('Positions (um)','FontSize', 14)
    ylabel('Counts','FontSize', 14)
    legend('FontSize', 14)
    title(['Position distribution: ' TypeOfPositionDistribution],'FontSize', 14)
    
    subplot(1,2,2)
    histogram(initialVelocities*1e3,NumberOfBins, 'LineStyle', 'none', 'DisplayName','Sampled')
    hold on
    plot(Solver.velocityGrid.*1e3, VelocityPDFForDisp, 'Color', [0.6350, 0.0780, 0.1840], 'LineStyle', '--', 'LineWidth', 1.5,'DisplayName','Predicted')
    xlabel('Velocities (mm/s)','FontSize', 14)
    ylabel('Counts','FontSize', 14)
    legend('FontSize', 14)
    title(['Velocity distribution: ' TypeOfVelocityDistribution ' with T = ' num2str(round(initialTemperature*1e9*1E1)/1E1) ' nK'], 'FontSize', 14)
    
    sgtitle('Sampling from a Distribution','FontSize', 18)
end