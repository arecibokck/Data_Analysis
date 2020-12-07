function ret = MaxwellBoltzmannDistribution(Temperature, Velocity)
    % ret = MaxwellBoltzmannDistribution(Temperature, Velocity)
    % computes the Probability density for a Cs133-atom to have a velocity
    % [m/s] at given Temperature [K]
    
    PhysicsConstants;
    m = Cs133Mass;
    k = BoltzmannConstant;
    %ret = 4*pi*sqrt((m/(2*pi*k*Temperature))^3) ...
    %     * Velocity.^2  .* exp((-m*Velocity.^2)/(2*k*Temperature));
    ret = sqrt((m/(2*pi*k*Temperature))^(1/2)) ...
         .* exp((-m*Velocity.^2)/(2*k*Temperature));
    ret = ret./sum(ret);
end
