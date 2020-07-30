clearvars
clc
poolobj = gcp('nocreate'); % Check if pool is open
if isempty(poolobj)
    parpool;
end
%% initializing the potential
PhysicsConstants;
Wavelength = 1064e-9;                               % VDT wavelength
Frequency  = 2*pi*SpeedOfLight/Wavelength;
detuningD1 = 2*pi*SpeedOfLight*(1/CsD1lambda-1/Wavelength);
detuningD2 = 2*pi*SpeedOfLight*(1/CsD2lambda-1/Wavelength);
% Trap parameters
Trap.k = 2*pi./Wavelength;
Trap.w0   = 50e-6;                                  % Beam waist
Trap.z_R  = pi.*(Trap.w0.^2)/Wavelength;            % Rayleigh length
Trap.w = @(z) Trap.w0.*sqrt(1+((z./Trap.z_R).^2));  % Axial waist
Trap.P = 500e-3;                                    % Beam power
fOscD1 = 0.344;                                     % D1 Absorption oscillator strength
fOscD2 = 0.714;                                     % Absorption oscillator strength
potentialContributionD1 = (fOscD1*CsD1Gamma/(2*pi*SpeedOfLight/CsD1lambda)^3)*(1/detuningD1+(1/(detuningD1+2*Frequency)));
potentialContributionD2 = (fOscD2*CsD2Gamma/(2*pi*SpeedOfLight/CsD2lambda)^3)*(1/detuningD2+(1/(detuningD2+2*Frequency)));
I0 = 2.*Trap.P./(pi.*(Trap.w0).^2);  % Single peak beam intensity, factor of 2 because two counter propagating beams
Imax = 2*I0; % Multiply another factor of 2 because of interference
Trap.U0 = -(3*pi*SpeedOfLight^2*Imax/2)*(potentialContributionD1+potentialContributionD2);
Trap.U0InTemperature = Trap.U0/BoltzmannConstant;
Trap.U0InFreq = Trap.U0/PlanckConstant;
% Plot the potential
x = -200*Wavelength:1e-2*Wavelength:200*Wavelength;
figure(1)
clf
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
set(gcf, 'defaultAxesColorOrder', [[0 0 0];[0 0 0]]);
Potential = @(x,Trap) Trap.U0 .* exp(-2 *((x/Trap.w0).^2));
%plot(x*1e6,-cumsum(force(x,Trap)), 'Color', colours{2}, 'LineWidth', 2) % plot the potential
plot(x*1e6, (Potential(x,Trap)./BoltzmannConstant) .*1e6, 'Color', colours{2}, 'LineWidth', 2) % plot the potential
xlabel('Position (um)','FontSize', 14)
ylabel('Potential (uK)','FontSize', 14)
xlim([min(x) max(x)]*1e6)
yyaxis right
plot(x*1e6,force(x,Trap)*1e23, 'Color', colours{4}, 'LineWidth', 2) % plot the potential
ylabel('Force (x 10^{-23} N)','FontSize', 14)
legend({'Potential', 'Force'},'FontSize', 14)
sgtitle(['Vertical Dipole Trap of depth = ' num2str(abs(Trap.U0InTemperature*1e6),'%.2f') ' uK'])
%% Simulation of trajectory of an atom allowed to oscillate in the trap
tRes        = 0.05e-6;              % Resolution for the ODE solver (s)
t0          = 0;                    % Starting time (s)
tf          = 2e-3;                % Final time (s)
tNumPoints  = floor(tf/tRes)+1;     % Number of sample points in time between t0 and tf
tspan = linspace(t0,tf,tNumPoints); % Solver calculates atom position for each of these timesteps in this time array
initialPositions = (1:1:10).*1e-6;
%initialPositions = (1:10:50).*1e-6;
%initialPositions = union((1:1:10), (10:10:100)).*1e-6;
NumberOfAtoms = length(initialPositions);
Xres = zeros(length(tspan),NumberOfAtoms);
Vres = zeros(length(tspan),NumberOfAtoms);  
progressbar  = parforNotifications();
progressbar.PB_start(NumberOfAtoms,'Message',['Computing evolution for ' num2str(NumberOfAtoms,'%.0f') ' atoms:']);
parfor Index = 1:NumberOfAtoms
    InitialConditions = [initialPositions(Index) 0]; % [initial position, initial velocity]
    [res] = ode5(@(t, x) odefcn(t, x, Trap), tspan, InitialConditions);
    Xres(:,Index) = res(:,1);
    Vres(:,Index) = res(:,2);
    progressbar.PB_iterate();
end
clear Index
%% Plotting
% Plot all trajectories 
figure(2)
clf
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
TimeForFullCompression = zeros(NumberOfAtoms,1);  
for Index = 1:NumberOfAtoms
    plot(tspan*1e3, Xres(:,Index).*1e6, 'HandleVisibility', 'Off');
    hold on
    ZC = ZeroX(tspan*1e3,Xres(:,Index).*1e6);
    TimeForFullCompression(Index) = ZC(1);    
end
MeanTime = mean(TimeForFullCompression);
line([min(TimeForFullCompression) min(TimeForFullCompression)],[-10 10],'Color',colours{1},'LineStyle','--')
line([MeanTime MeanTime],[-10 10],'Color',colours{1},'LineWidth',1.5)
line([max(TimeForFullCompression) max(TimeForFullCompression)],[-10 10],'Color',colours{1},'LineStyle','--')
clear Index
sgtitle(['Trajectory of ' num2str(NumberOfAtoms) ' atoms at different starting positions in the VDT of depth: ' num2str(abs(Trap.U0InTemperature*1e6),'%.2f') ' uK']);
ylabel('Position (um)','FontSize', 14)
xlabel('Time (ms)','FontSize', 14)
legend({['Minimum quarter period (' num2str(min(TimeForFullCompression),'%.3f') ' ms)'], ['Mean quarter period (' num2str(MeanTime,'%.3f') ' ms)'], ['Maximum quarter period (' num2str(max(TimeForFullCompression),'%.3f') ' ms)']}, 'FontSize', 14)
grid on
hold off
figure(3)
clf
set(gcf, 'defaultAxesColorOrder', [colours{1};colours{5}]);
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
yyaxis left
plot(initialPositions*1e6, TimeForFullCompression, '--o', 'Color',colours{1})
sgtitle(['For atoms starting at different positions in the VDT of depth: ' num2str(abs(Trap.U0InTemperature*1e6),'%.2f') ' uK']);
xlabel('Starting Position (um)','FontSize', 14)
ylabel('Quarter period [ms]','FontSize', 14)
grid on
yyaxis right
plot(initialPositions*1e6, 1E3./(4.*TimeForFullCompression), '--o', 'Color',colours{5})
ylabel('Trapping Frequency [Hz]','FontSize', 14)
legend({'Quarter period','Trapping Frequency'}, 'FontSize', 14)
%% Save
prompt = 'Save trajectories? Enter "true" or "false": ';
SaveTrajectories = input(prompt);
if isempty(SaveTrajectories)
    SaveTrajectories = false;
end
if SaveTrajectories
    save('Trajectories.mat','tspan','Xres')
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