function plotEnvelopes(Solver, UptoInitPos)
    %%
    f_h = DQSIMhelper.getFigureByTag('EnvelopesAndOptimalWaitTimes');
    set(groot,'CurrentFigure',f_h);
    a_h = get(f_h, 'CurrentAxes');
    if ~isempty(get(a_h, 'Children')) 
        clf(f_h);
    end
    f_h.Name = 'EnvelopesAndOptimalWaitTimes';
    f_h.Units = 'pixels';
    
    f_h.Position = [Solver.FigurePosition 760 600];
    
    positions  = Solver.simulationResults(:,:,1);
    
    tspan = Solver.timeSpan;
    
    initialPositions = positions(1,:);
    [initialPositions, sortIdx] = sort(initialPositions, 'ascend');
    positions = positions(:,sortIdx);
    envelopes = zeros(size(positions));
    [~,OriginIndex]=min(abs(round(initialPositions.*1e6)));
    for idx = 1:size(positions,2)-1
        FirstPosition = round(positions(1,idx).*1e6);
        if  FirstPosition > 0
            Limit = find(round(initialPositions.*1e6)==FirstPosition, 1, 'first');
            envelopes(:,idx) = max(abs(positions(:,OriginIndex:Limit)),[],2).*1e6;
        elseif FirstPosition < 0
            Limit = find(round(initialPositions.*1e6)==FirstPosition, 1, 'first');
            envelopes(:,idx) = max(abs(positions(:,Limit:OriginIndex)),[],2).*1e6;
        end
    end
    idx = find(round(initialPositions.*1e6)==UptoInitPos, 1, 'first');
    idx = idx(1);
    plot(tspan*1e3, envelopes(:,idx), 'LineWidth', 5, 'Color',[0.8500, 0.3250, 0.0980]);
    hold on
    plot(tspan*1e3, abs(positions).*1e6,'HandleVisibility', 'Off');
    xlabel('Time (ms)','FontSize', 14)
    ylabel('Position (um)','FontSize', 14)
    legend({['Upto ' num2str(initialPositions(idx)*1e6) ' um']}, 'FontSize', 14)
    sgtitle('Optimal waiting times for different capture ranges');
end