clearvars
clc
poolobj = gcp('nocreate'); % Check if pool is open
if isempty(poolobj)
    parpool;
end
global Power c potentialContributionD1 potentialContributionD2 PositionArray tRes t0 tf;
PhysicsConstants;
Wavelength = 1064e-9;                               % VDT wavelength
c          = SpeedOfLight;
Frequency  = 2*pi*c/Wavelength;
detuningD1 = 2*pi*c*(1/CsD1lambda-1/Wavelength);
detuningD2 = 2*pi*c*(1/CsD2lambda-1/Wavelength);
% Trap parameters
waist_guess = 50e-6;                                % Beam waist
Power      = 500e-3;                                % Beam power
fOscD1     = 0.344;                                 % D1 Absorption oscillator strength
fOscD2     = 0.714;                                 % Absorption oscillator strength
potentialContributionD1 = (fOscD1*CsD1Gamma/(2*pi*c/CsD1lambda)^3)*(1/detuningD1+(1/(detuningD1+2*Frequency)));
potentialContributionD2 = (fOscD2*CsD2Gamma/(2*pi*c/CsD2lambda)^3)*(1/detuningD2+(1/(detuningD2+2*Frequency)));
PositionArray = 1:10:100;
tRes        = 0.5e-6;                              % Resolution for the ODE solver (s)
t0          = 0;                                    % Starting time (s)
tf          = 20e-3;                                 % Final time (s)

fcn = @(w) sum(diff(Optimizer(w)));
lb = 50e-6;
ub = 150e-6;
optimumwaist = fmincon(fcn, waist_guess,[],[],[],[],lb,ub); 
% Local minimum was found that satisfies the constraints.
% This minimum was 64.261e-06 m. <-----Took a long time!
% Optimization completed because the objective function is non-decreasing in 
% feasible directions, to within the value of the optimality tolerance,
% and constraints are satisfied to within the value of the constraint tolerance.

function TimesForFullCompression = Solver(tRes, t0, tf, Trap)
    global PositionArray
    tNumPoints  = floor(tf/tRes)+1;     % Number of sample points in time between t0 and tf
    tspan = linspace(t0,tf,tNumPoints); % Solver calculates atom position for each of these timesteps in this time array
    initialPositions = PositionArray.*1e-6;
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
    TimesForFullCompression = zeros(NumberOfAtoms,1);  
    for Index = 1:NumberOfAtoms
        ZC = ZeroX(tspan*1e3,Xres(:,Index).*1e6);
        TimesForFullCompression(Index) = ZC(1); 
    end
end
function ret = Optimizer(Waist)
    global Power c potentialContributionD1 potentialContributionD2 tRes t0 tf;
    % Trap parameters
    Trap.w0   = Waist;                                  % Beam waist
    Trap.P    = Power;                                  % Beam power
    I0 = 2.*Trap.P./(pi.*(Trap.w0).^2);                 % Single peak beam intensity, factor of 2 because two counter propagating beams
    Imax = 2*I0;                                        % Multiply another factor of 2 because of interference
    Trap.U0 = -(3*pi*c^2*Imax/2)*(potentialContributionD1+potentialContributionD2);
    ret = Solver(tRes, t0, tf, Trap);
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