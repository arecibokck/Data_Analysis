function plotPotential_private(this, Type)
    
    temp = this.potentialType;
    this.potentialType = Type;
    
    f_h = DQSIMhelper.getFigureByTag('potentials');
    figure(f_h.Number);
    f_h.Name = ['Potentials ' this.compressionDirection ' compression'];
    f_h.Units = 'pixel';
    
    f_h.Position = [this.FigurePosition 760 600];
    
    yyaxis left
    U = (this.potential(this.positionGrid)./this.physicalConstants.BoltzmannConstant) .*1e6;
    p_1 = plot(this.positionGrid*1e6, U, 'Color', [0, 0.4470, 0.7410], 'LineStyle', '-', 'LineWidth', 2); % plot the potential
    xlabel([CapitalizeString(this.compressionDirection) ' position (um)'],'FontSize', 14)
    ylabel('Potential (uK)','FontSize', 14)
    xlim([min(this.positionGrid) max(this.positionGrid)]*1e6) 
    
    yyaxis right
    F = this.force(this.positionGrid)*1e23;
    p_2 = plot(this.positionGrid*1e6, F, 'Color', [0.8500, 0.3250, 0.0980], 'LineStyle', '--', 'LineWidth', 2); % plot the force
    ylabel('Force (x 10^{-23} N)','FontSize', 14)
    
    legend([p_1 p_2], {'Potential', 'Force'},'FontSize', 14)
    sgtitle(['Dipole Trap of depth = ' num2str(this.potentialDepthInTemp,'%.2f') ' uK'])
    
    this.potentialType = temp;
end

function str = CapitalizeString(str)
    str(1) = upper(str(1));
end
