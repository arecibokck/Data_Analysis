function plotConvolvedDistributionSampling(Solver)
    % - plot Distribution
    f_h = DQSIMhelper.getFigureByTag('ConvolvedDistribution');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children'))
        clf(f_h);
    end
    f_h.Name = 'Sampling from a convolved distribution';
    f_h.Units = 'pixels';

    f_h.Position = [Solver.FigurePosition 1.4293e+03 824];
    pc = Solver.physicalConstants;

    n = Solver.numberOfAtoms;
    initialTemperature = Solver.InitialDistributionParameters.initialTemperature;
    velocities = Solver.velocityGrid;
    
    TrapDepthAtRelease = sqrt(Solver.InitialDistributionParameters.FractionOfInitialPotential) *  LatticeProperties.estimateTrapDepthFromHeatingSidebandFreqDT1DT3(Solver.trapFreqHDTy);
    results = LatticeProperties.estimateTrapFreq33FromTrapDepth(0,0,0,0,TrapDepthAtRelease);
    TrapFrequencyAtRelease = results(1);
    DeltaVQuantumMechanical = sqrt(0.5*pc.Cs133Mass*pc.PlanckConstantReduced*TrapFrequencyAtRelease)./pc.Cs133Mass;
    ProbDFArray = cell(1,3);
    ProbDFArray{1}    = CompressionHelper.MaxwellBoltzmannDistribution(initialTemperature, velocities);
    ProbDFArray{2}    = CompressionHelper.GaussianDistribution(0, DeltaVQuantumMechanical, velocities);
    ProbDFArray{3}    = conv(ProbDFArray{1}, ProbDFArray{2});
    ProbDFName = {'Maxwell-Boltzmann','Gaussian','Convolution'};
    for kk = 1:3
        ProbDF = ProbDFArray{kk};
        ProbDF     = ProbDF./sum(ProbDF);
        subplot(1,3,kk);
        switch kk
            case 1
                stairs(velocities.*1e3, ProbDF, 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 2)
                xlabel('Velocities(mm/s)','FontSize', 14)
                ylabel('Probability','FontSize', 14)
                yyaxis right
                plot(velocities.*1e3, cumsum(ProbDF))
            case 2
                MaxVelocity          = max(velocities);
                MinVelocity          = min(velocities);
                RangeExtensionFactor = 0.25;
                velocityVals = linspace(MinVelocity - RangeExtensionFactor, MaxVelocity + RangeExtensionFactor, length(velocities));
                ProbDF       = CompressionHelper.GaussianDistribution(0, DeltaVQuantumMechanical, velocityVals);
                ProbDF       = ProbDF./sum(ProbDF);
                stairs(velocityVals.*1e3, ProbDF, 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 2)
                xlabel('Velocities(mm/s)','FontSize', 14)
                ylabel('Probability','FontSize', 14)
                yyaxis right
                plot(velocityVals.*1e3, cumsum(ProbDF))
            case 3
                MaxVelocity = max(velocities);
                MinVelocity = min(velocities);
                NumberOfBins = 10;
                velocityVals = MinVelocity: (MaxVelocity-MinVelocity)/(length(ProbDF)-1):MaxVelocity;
                VelocityPDFForDisp = ProbDF * n * numel(ProbDF)/NumberOfBins;
                initialVelocities = CompressionHelper.drawSamplesFromDistribution(n, velocityVals, ProbDF);
                histogram(initialVelocities.*1e3, NumberOfBins, 'FaceAlpha', 0.1)
                hold on
                plot(velocityVals.*1e3, VelocityPDFForDisp, 'Color', [0.6350, 0.0780, 0.1840], 'LineStyle', '--', 'LineWidth', 1.5,'DisplayName','Predicted')
                xlabel('Velocities(mm/s)','FontSize', 14)
                ylabel('Frequency','FontSize', 14)
                legend({'Sampled Velocities', 'PDF'})
        end
        title(ProbDFName{kk},'FontSize', 14)
    end
    sgtitle('Sampling from a convolution of distributions','FontSize', 18)
end
