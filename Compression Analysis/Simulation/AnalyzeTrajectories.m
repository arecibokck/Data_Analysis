clearvars
clc
poolobj = gcp('nocreate'); % Check if pool is open
if isempty(poolobj)
    parpool;
end
load('Trajectories_N1000.mat')
NumberOfAtoms = size(Xres, 2);
colours     = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
tRes        = 0.5e-6;               % Resolution for the ODE solver (s)
t0          = 0;                    % Starting time (s)
tf          = 10e-3;                % Final time (s)
tNumPoints  = floor(tf/tRes)+1;     % Number of sample points in time between t0 and tf
%% - Plotting
% figure(1)
% clf
% set(gcf, 'Units', 'normalized');
% set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
% c = linspace(0,1,NumberOfAtoms);
% F = struct('cdata',[],'colormap',[]);
% frame = 0;
% for Time = 1:50:tNumPoints
%     positions  = Xres(Time,:).*1e6;
%     velocities = Vres(Time,:).*1e3;
%     PhaseSpaceMatrix = [positions(:), velocities(:)];
%     hist2d(positions(:),  velocities(:), ...
%         'nbins', 20, ...
%         'PositionLimits', [-round(max(max(Xres))*1e6, 1) round(max(max(Xres))*1e6, 1)],...
%         'VelocityLimits', [-round(max(max(Vres))*1e3, 1) round(max(max(Vres))*1e3, 1)],...
%         'CountDensity', false);
%     colorbar
%     %         points = scatter(positions, velocities, 20, abs(velocities),'filled', 'HandleVisibility', 'Off');
%     %         colormap(jet(size(velocities,2)))
%     xlabel('Position (\mum)','FontSize', 14)
%     ylabel('Velocity (mm/s)','FontSize', 14)
%     xlim([-round(max(max(Xres))*1e6, 1) round(max(max(Xres))*1e6, 1)])
%     ylim([-round(max(max(Vres))*1e3, 1) round(max(max(Vres))*1e3, 1)])
%     %sgtitle(['Ground State pop = ' num2str(GroundStatePopulation) ' --> (Uniform) Initial Temperature = ' num2str(initialTemperature*1e9, '%.1f') ' nK']);
%     grid on
%     hold on
%     frame = frame+1;
%     F (frame) = getframe (gcf);
%     drawnow
%     %         if Time~=tNumPoints
%     %             points.MarkerFaceAlpha = .05;
%     %         else
%     %             points.MarkerFaceAlpha = 1;
%     %         end
% end
% 
% clear Index
% hold off
% writerObj = VideoWriter (['PhaseEvolution_N' num2str(NumberOfAtoms) '.avi']);
% writerObj.FrameRate = 15;
% writerObj.Quality = 100;
% open (writerObj);
% for i = 1: length (F)
%     frame = F (i);
%     writeVideo (writerObj, frame);
% end
% % close the writer object
% close (writerObj);
%%
figure(2)
clf
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
TimeForFullCompression = zeros(NumberOfAtoms,1);  
for Index = 1:NumberOfAtoms
    plot(tspan*1e3, Xres(:,Index).*1e6, 'HandleVisibility', 'Off');
    hold on
    ZC = findAllZeroCrossings(tspan*1e3,Xres(:,Index).*1e6);
    TimeForFullCompression(Index) = ZC(1);    
end
clear Index
MeanTime = mean(TimeForFullCompression);
line([min(TimeForFullCompression) min(TimeForFullCompression)],[-100 100],'Color',colours{1},'LineStyle','--')
%line([MeanTime MeanTime],[-100 100],'Color',colours{1},'LineWidth',1.5)
line([max(TimeForFullCompression) max(TimeForFullCompression)],[-100 100],'Color',colours{1},'LineStyle','--')
sgtitle(['Trajectory of ' num2str(NumberOfAtoms) ' atoms at different starting positions in the VDT']);
ylabel('Position (um)','FontSize', 14)
xlabel('Time (ms)','FontSize', 14)
%legend({['Minimum quarter period (' num2str(min(TimeForFullCompression),'%.3f') ' ms)'], ['Mean quarter period (' num2str(MeanTime,'%.3f') ' ms)'], ['Maximum quarter period (' num2str(max(TimeForFullCompression),'%.3f') ' ms)']}, 'FontSize', 14)
legend({['Minimum quarter period (' num2str(min(TimeForFullCompression),'%.3f') ' ms)'], ['Maximum quarter period (' num2str(max(TimeForFullCompression),'%.3f') ' ms)']}, 'FontSize', 14)
grid on
hold off
%%
initialPositions = Xres(1,:);
figure(3)
clf
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
set(gcf, 'defaultAxesColorOrder', [colours{1};colours{5}]);
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
%yyaxis left
[initialPositions, sortIdx] = sort(initialPositions, 'ascend');
TimeForFullCompression = TimeForFullCompression(sortIdx);
plot(initialPositions*1e6, TimeForFullCompression, '--o', 'Color',colours{1})
%sgtitle(['For atoms starting at different positions in the VDT of depth: ' num2str(abs(Trap.U0InTemperature*1e6),'%.2f') ' uK']);
sgtitle('For atoms starting at different positions in the VDT');
xlabel('Starting Position (um)','FontSize', 14)
ylabel('Quarter period [ms]','FontSize', 14)
grid on
% yyaxis right
% plot(initialPositions*1e6, 1E3./(4.*TimeForFullCompression), '--o', 'Color',colours{5})
% ylabel('Trapping Frequency [Hz]','FontSize', 14)
% legend({'Quarter period','Trapping Frequency'}, 'FontSize', 14)
%%
figure(4);
clf;
envelopes = zeros(size(Xres));
optimalwaitingtimes  = zeros(size(Xres,2)-1,1);
for idx = 1:size(Xres,2)-1
    envelopes(:,idx) = max(abs(Xres(:,1:idx+1)),[],2).*1e6;
    [~, index] = min(envelopes(:,idx));
    optimalwaitingtimes(idx) = tspan(index)*1e3;
end
InitPos = 81; %in um
idx = find(round(initialPositions.*1e6)==InitPos);
plot(tspan*1e3, envelopes(:,idx), 'LineWidth', 5, 'Color',colours{2});
hold on
plot(tspan*1e3, abs(Xres).*1e6,'HandleVisibility', 'Off');
xlabel('Time (ms)','FontSize', 14)
ylabel('Position (um)','FontSize', 14)
legend({['Upto ' num2str(initialPositions(idx)*1e6) ' um: \tau_{wait}^{opt} ~ ' num2str(optimalwaitingtimes(12), '%.2f') ' ms']}, 'FontSize', 14)
sgtitle('Optimal waiting times for different capture ranges');
%%
figure(5);
clf
plot(initialPositions(2:end-1)*1e6,optimalwaitingtimes(1:size(Xres, 2)-2), '--o', 'LineWidth', 2)
xlabel('Capture Range (um)','FontSize', 14)
ylabel('Optimal Waiting Time (ms)','FontSize', 14)
sgtitle('Optimal waiting times for different capture ranges');
function ZC  = findAllZeroCrossings(x,y)
% ZC = findAllZeroCrossings(x,y) finds all Zero-crossing of the function y = f(x)
%
% 
% Remark:
%   findAllZeroCrossings has been modified so as to avoid the following error "Error using griddedInterpolant
%   The grid vectors are not strictly monotonic increasing."
%
%   The following modifications are based off of this MatLab forum answer:
%   https://www.mathworks.com/matlabcentral/answers/283718-grid-vectors-not-strictly-monotonic-increasing
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns Approximate Zero-Crossing Indices Of Argument Vector
zxidx = zci(y);
if ~isempty(zxidx)
    for k1 = 1:numel(zxidx)
        idxrng = max([1 zxidx(k1)-1]):min([zxidx(k1)+1 numel(y)]);
        xrng = x(idxrng);
        yrng = y(idxrng);
        % Beginning of findAllZeroCrossings2 modifications. The naming conventions follow
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
function hist2d(x,y,varargin)
    
    p = inputParser;
    addRequired(p,  'Positions',            @isnumeric)
    addRequired(p,  'Velocities',           @isnumeric)
    addParameter(p, 'nbins',           10,  @isscalar)
    addParameter(p, 'PositionLimits',  [-60 60],  @(x) isDataArray(x) && isnumeric(x))
    addParameter(p, 'VelocityLimits',  [-80 80],  @(x) isDataArray(x) && isnumeric(x))
    addParameter(p, 'CountDensity', false,  @islogical)
    parse(p,x,y,varargin{:})
    
    x = p.Results.Positions;
    y = p.Results.Velocities;
    nbins = p.Results.nbins;
    plim  = p.Results.PositionLimits;
    vlim  = p.Results.VelocityLimits;
    cdflag  = p.Results.CountDensity;
    
    xedges = linspace(min(plim),max(plim),nbins+1);
    yedges = linspace(min(vlim),max(vlim),nbins+1);
    
    % computes bins vectors (as middle of each edges couples)
    xbins = mean(cat(1,xedges(1:end-1),xedges(2:end)));
    ybins = mean(cat(1,yedges(1:end-1),yedges(2:end)));

    % computes bins width vectors and area matrix
    xbw = diff(xedges);
    ybw = diff(yedges);
    [xx,yy] = meshgrid(xbw,ybw);
    a = xx.*yy;
    
    % initiate the result matrix
    n = zeros(length(ybins),length(xbins));
    % main loop to fill the matrix with element counts
    for i = 1:size(n,1)
        k = find(y >= yedges(i) & y < yedges(i+1));
        for j = 1:size(n,2)
            n(i,j) = length(find(x(k) >= xedges(j) & x(k) < xedges(j+1)));
        end
    end
    
    % normalize options
    if cdflag
        n = n./a; % count density
    end
    imagesc(xbins,ybins,n)
    %hold on
    %plot(x,y,'.k','MarkerSize',10)
    %hold off	
end