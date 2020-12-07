function runClassicalSimulation(this)
    % - perform the simulation
    
    tRes        = this.timeStep;        % Resolution for the ODE solver (s)
    t0          = 0;                    % Starting time (s)
    tf          = this.finalTime;       % Final time (s)
    tNumPoints  = floor(tf/tRes)+1;     % Number of sample points in time between t0 and tf
    this.timeSpan = linspace(t0,tf,tNumPoints); % Solver calculates atom position for each of these timesteps in this time array
    
    n = this.numberOfAtoms;
    
    initialPositions = this.initialPositions;
    initialVelocities = this.initialVelocities;
    
    positions = zeros(length(this.timeSpan),n);
    velocities = zeros(length(this.timeSpan),n);
    
    progressbar  = CompressionHelper.parforNotifications();
    progressbar.PB_start(n,'Message',['Computing evolution for ' num2str(n,'%.0f') ' atoms:']);
    
    parfor Index = 1:n
        if ~isscalar(initialVelocities)
            InitialConditions = [initialPositions(Index) initialVelocities(Index)]; % [initial position, initial velocity]
        else
            InitialConditions = [initialPositions(Index) initialVelocities];
        end
        [res] = CompressionHelper.ode5(@(t, x) this.odefcn(t, x), this.timeSpan, InitialConditions);
        positions(:,Index) = res(:,1);
        velocities(:,Index) = res(:,2);
        progressbar.PB_iterate();
    end
    clear Index
    
    this.simulationResults = cat(3, positions, velocities);
    
    % Save
    if this.DoSave
        save(['Trajectories_N' num2str(n) '.mat'],'tspan','positions')
        save(['Trajectories_N' num2str(n) '.mat'],'velocities', '-append')
    end
    
end