function plotPhaseSpaceEvolution(Solver, NumberOfBins)
    f_h = DQSIMhelper.getFigureByTag('PhaseSpaceEvolution');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'Phase-Space Evolution';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 990 770];
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
    F = struct('cdata',[],'colormap',[]);
    frame = 0;
    
    GroundStatePopulation = Solver.InitialDistributionParameters.GroundStatePopulation;
    initialTemperature = Solver.InitialDistributionParameters.initialTemperature;
    FinalTrapDepthInUnitsOfRecoilEnergy = Solver.InitialDistributionParameters.FinalTrapDepthInUnitsOfRecoilEnergy;
    
    tRes        = Solver.timeStep;        % Resolution for the ODE solver (s)
    tf          = Solver.finalTime;       % Final time (s)
    tNumPoints  = floor(tf/tRes)+1;     % Number of sample points in time between t0 and tf
    tspan = Solver.timeSpan;
    
    n = Solver.numberOfAtoms;
    pc = Solver.physicalConstants;
    allPositions  = Solver.simulationResults(:,:,1);
    allVelocities = Solver.simulationResults(:,:,2);
    
    EscapeVelocities = sqrt(2.*abs(Solver.potential(Solver.positionGrid))./pc.Cs133Mass)*1e3;
    
    for Time = 1:5:tNumPoints
        positions  = allPositions(Time,:).*1e6;
        velocities = allVelocities(Time,:).*1e3;
        sb1 = subplot(4,4,[1,9]);
        h1 = histogram(velocities(:), NumberOfBins,'Normalization','probability');
        h1.FaceAlpha =  0.1;
        h1.FaceColor = colours{2};
        ylabel('Velocity Distribution','FontSize', 12)
        [MaxValue,~] = max(h1.Values);
        set(sb1, 'Box', 'off', 'Color', 'none')
        set(sb1, 'xlim', [-round(max(allVelocities(:))*1e3, 1) round(max(allVelocities(:))*1e3, 1)])
        set(sb1, 'ylim', [0 MaxValue])
        set(sb1, 'XDir', 'reverse')
        camroll(sb1,90)
        hold off
        sb2 = subplot(4,4,[2.15,12]);
        Plotting.plot2DHistogram(positions(:),  velocities(:), ...
            'nbins', NumberOfBins, ...
            'PositionLimits', [-round(max(allPositions(:))*1e6, 1) round(max(allPositions(:))*1e6, 1)], ...
            'VelocityLimits', [-round(max(allVelocities(:))*1e3, 1) round(max(allVelocities(:))*1e3, 1)], ...
            'CountDensity', false);
        colorbar
        hold on
        plot(Solver.positionGrid*1e6, EscapeVelocities, 'Color', [1 1 1], 'LineStyle', '--')
        plot(Solver.positionGrid*1e6, -EscapeVelocities, 'Color', [1 1 1], 'LineStyle', '--')
        text(round(min(allPositions(:))*1e6, 1)+2, -round(max(allVelocities(:))*1e3, 1)+5, ['Time : ' num2str(tspan(Time)*1e3) ' ms'], 'Color', [1 1 1], 'FontSize', 14)
        xlabel('Position (\mum)','FontSize', 14)
        ylabel('Velocity (mm/s)','FontSize', 14)
        xlim([-round(max(allPositions(:))*1e6, 1) round(max(allPositions(:))*1e6, 1)])
        ylim([-round(max(allVelocities(:))*1e3, 1) round(max(allVelocities(:))*1e3, 1)])
        grid on
        sb3 = subplot(4,4,[14.15,15.7], 'color', 'none');
        h2 = histogram(positions(:), NumberOfBins,'Normalization','probability');
        h2.FaceAlpha =  0.1;
        h2.FaceColor = colours{2};
        ylabel('Position Distribution','FontSize', 12)
        [MaxValue,~] = max(h2.Values);
        set(sb3, 'Box', 'off', 'Color', 'none')
        set(sb3, 'xlim', [-round(max(allPositions(1,:))*1e6, 1) round(max(allPositions(1,:))*1e6, 1)])
        set(sb3, 'ylim', [0 MaxValue])
        hold off
        sgtitle(['Ground State Pop. = ' num2str(GroundStatePopulation) ' --> (Uniform) Initial Temperature = ' num2str(initialTemperature*1e6, '%.1f') ' \muK; Release Depth = ' num2str(FinalTrapDepthInUnitsOfRecoilEnergy) ' E_R']);
        frame = frame+1;
        F (frame) = getframe (gcf);
        drawnow
    end
    
    clear Index
    hold off
    writerObj = VideoWriter (['PhaseEvolution_N' num2str(n) '.avi']);
    writerObj.FrameRate = 10;
    writerObj.Quality = 100;
    open (writerObj);
    for i = 1: length (F)
        frame = F (i);
        writeVideo (writerObj, frame);
    end
    % close the writer object
    close (writerObj);
end
