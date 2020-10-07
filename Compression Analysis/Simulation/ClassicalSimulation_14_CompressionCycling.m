clearvars
clc
poolobj = gcp('nocreate'); % Check if pool is open
if isempty(poolobj)
    parpool;
end
PhysicsConstants;
Wavelength = 1064e-9;                               % VDT wavelength
c          = SpeedOfLight;
Frequency  = 2*pi*c/Wavelength;
detuningD1 = 2*pi*c*(1/CsD1lambda-1/Wavelength);
detuningD2 = 2*pi*c*(1/CsD2lambda-1/Wavelength);
% Trap parameters
waist_guess = 50e-6;                                % Beam waist
Trap.P     = 500e-3;                                % Beam power
fOscD1     = 0.344;                                 % D1 Absorption oscillator strength
fOscD2     = 0.714;                                 % Absorption oscillator strength
potentialContributionD1 = (fOscD1*CsD1Gamma/(2*pi*c/CsD1lambda)^3)*(1/detuningD1+(1/(detuningD1+2*Frequency)));
potentialContributionD2 = (fOscD2*CsD2Gamma/(2*pi*c/CsD2lambda)^3)*(1/detuningD2+(1/(detuningD2+2*Frequency)));
AtomPositionDistribution = (1:2:80).*1e-6;
spread = max(AtomPositionDistribution).*1e6;
%AtomVelocityDistribution = zeros(1,length(AtomPositionDistribution));
AtomVelocityDistribution = (500:41:2130).*1e-6;
tRes        = 0.5e-6;                               % Resolution for the ODE solver (s)
t0          = 0;                                    % Starting time (s)
WaitTime    = 2e-3;                                 % Final time (s)
CompressionCycles = 3;                              % CompressionCycles
NumberOfAtoms = length(AtomPositionDistribution);
tNumPoints  = floor(CompressionCycles*WaitTime/tRes)+CompressionCycles;             % Number of sample points in time between t0 and final time
tspan = linspace(t0,CompressionCycles*WaitTime,tNumPoints);
Trap.w0     = 50e-6;                                % Beam waist
I0 = 2.*Trap.P./(pi.*(Trap.w0).^2);                 % Single peak beam intensity, factor of 2 because two counter propagating beams
Imax = 2*I0;                                        % Multiply another factor of 2 because of interference
Trap.U0 = -(3*pi*c^2*Imax/2)*(potentialContributionD1+potentialContributionD2);
Trajectories = zeros(tNumPoints, NumberOfAtoms);
for idx = 1:CompressionCycles
    [CurrentPositionDistribution,CurrentVelocityDistribution] = Solver(tRes, t0, WaitTime, Trap, AtomPositionDistribution, AtomVelocityDistribution);
    Trajectories((idx-1)*size(CurrentPositionDistribution,1)+1:(idx*size(CurrentPositionDistribution,1)),:) = CurrentPositionDistribution;  
    AtomPositionDistribution = CurrentPositionDistribution(end,:);
    %AtomVelocityDistribution = CurrentVelocityDistribution(end,:);
end
%% Plotting
% Plot all trajectories 
figure(1)
clf
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
for Index = 1:NumberOfAtoms
    plot(tspan*1e3, Trajectories(:,Index).*1e6, 'HandleVisibility', 'Off');
    hold on
end
clear Index
for Index = 1:CompressionCycles
    line([Index*(WaitTime*1e3) Index*(WaitTime*1e3)],[-100 100],'Color',colours{1},'LineStyle','--')
end
clear Index
sgtitle(['Trajectory of ' num2str(NumberOfAtoms) ' atoms up to ' num2str(spread, '%.f') ' um for ' num2str(CompressionCycles) ' cycles']);
ylabel('Position (um)','FontSize', 14)
xlabel('Time (ms)','FontSize', 14)
grid on
hold off
function [Xres, Vres] = Solver(tRes, t0, tf, Trap, AtomPositionDistribution, AtomVelocityDistribution)
    tNumPoints  = floor(tf/tRes)+1;     % Number of sample points in time between t0 and tf
    tspan = linspace(t0,tf,tNumPoints); % Solver calculates atom position for each of these timesteps in this time array
    NumberOfAtoms = length(AtomPositionDistribution);
    Xres = zeros(length(tspan),NumberOfAtoms);
    Vres = zeros(length(tspan),NumberOfAtoms);  
    progressbar  = parforNotifications();
    progressbar.PB_start(NumberOfAtoms,'Message',['Computing evolution for ' num2str(NumberOfAtoms,'%.0f') ' atoms:']);
    parfor Index = 1:NumberOfAtoms
        InitialConditions = [AtomPositionDistribution(Index) AtomVelocityDistribution(Index)]; % [initial position, initial velocity]
        [res] = ode5(@(t, x) odefcn(t, x, Trap), tspan, InitialConditions);
        Xres(:,Index) = res(:,1);
        Vres(:,Index) = res(:,2);
        progressbar.PB_iterate();
    end
    clear Index
end
%% helper functions
function ret = odefcn(~, x, Trap)
    PhysicsConstants;
    ret = zeros(2,1);
    ret(1) = x(2); % Velocity
    ret(2) = force(x(1), Trap)/Cs133Mass;
end % - ODE to solve
function ret = force(x, Trap)
    ret = Trap.U0 * 4 * x/(Trap.w0.^2) .* exp(-2 *((x/Trap.w0).^2));
end % - Dipole force given a Gaussian potential  
function ZC = ZeroX(x,y)
% ZeroX has been modified so as to avoid the following error "Error using griddedInterpolant
% The grid vectors are not strictly monotonic increasing."
%
% The following modifications are based off of this MatLab forum answer:
% https://www.mathworks.com/matlabcentral/answers/283718-grid-vectors-not-strictly-monotonic-increasing
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns Approximate Zero-Crossing Indices Of Argument Vector
zxidx = zci(y);
if ~isempty(zxidx)
    for k1 = 1:numel(zxidx)
        idxrng = max([1 zxidx(k1)-1]):min([zxidx(k1)+1 numel(y)]);
        xrng = x(idxrng);
        yrng = y(idxrng);
        % Beginning of ZeroX2 modifications. The naming conventions follow
        % those in the referenced MatLab forum, except that "X" is "yrng" and
        % "Y" is "xrng".
        [yrng2, ~, jyrng] = unique(yrng); %yrng is a new array containing the unique values of yrng. jyrng contains the indices in yrng that correspond to the original vector. yrng = yrng2(jyrng)
        xrng2 = accumarray(jyrng, xrng, [], @mean); %This function creates a new array "xrng2" by applying the function "@mean" to all elements in "xrng" that have identical indices in "jyrng". Any elements with identical X values will have identical indices in jyrng. Thus, this function creates a new array by averaging values with identical X values in the original array.
        ZC(k1) = interp1( yrng2(:), xrng2(:), 0, 'linear', 'extrap' );
    end
else
    warning('No zero crossings found!')
    ZC = nan;
end
end