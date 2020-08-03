load('Trajectories.mat')
NumberOfAtoms = size(Xres, 2);
initialPositions = Xres(1,1:size(Xres, 2));
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};

%% - Plotting
figure(1)
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
figure(2)
clf
%set(gcf, 'defaultAxesColorOrder', [colours{1};colours{5}]);
set(gcf, 'Units', 'normalized');
set(gcf, 'OuterPosition', [0.5 0.5 0.25 0.45]);
%yyaxis left
plot(initialPositions*1e6, TimeForFullCompression, '--o', 'LineWidth', 2, 'Color',colours{1})
sgtitle('For atoms starting at different positions in the VDT');
xlabel('Starting Position (um)','FontSize', 14)
ylabel('Quarter period [ms]','FontSize', 14)
grid on
%yyaxis right
%plot(initialPositions*1e6, 1E3./(4.*TimeForFullCompression), '--o', 'Color',colours{5})
%ylabel('Trapping Frequency [Hz]','FontSize', 14)
%legend({'Quarter period','Trapping Frequency'}, 'FontSize', 14)
figure(3);
clf;
envelopes = zeros(size(Xres));
optimalwaitingtimes  = zeros(size(Xres,2)-1,1);
for idx = 1:size(Xres,2)-1
    envelopes(:,idx) = max(abs(Xres(:,1:idx+1)),[],2).*1e6;
    [~, index] = min(envelopes(:,idx));
    optimalwaitingtimes(idx) = tspan(index)*1e3;
end
InitPos = 80; %in um
idx = find(round(initialPositions.*1e6)==InitPos);
plot(tspan*1e3, envelopes(:,idx), 'LineWidth', 5, 'Color',colours{2});
hold on
plot(tspan*1e3, abs(Xres).*1e6,'HandleVisibility', 'Off');
xlabel('Time (ms)','FontSize', 14)
ylabel('Position (um)','FontSize', 14)
legend({['Upto ' num2str(initialPositions(idx)*1e6) ' um: \tau_{wait}^{opt} ~ ' num2str(optimalwaitingtimes(12), '%.2f') ' ms']}, 'FontSize', 14)
sgtitle('Optimal waiting times for different capture ranges');
figure(4);
clf
plot(initialPositions(2:end-1)*1e6,optimalwaitingtimes(1:size(Xres, 2)-2), '--o', 'LineWidth', 2)
xlabel('Capture Range (um)','FontSize', 14)
ylabel('Optimal Waiting Time (ms)','FontSize', 14)
sgtitle('Optimal waiting times for different capture ranges');
function ZC = ZeroX(x,y)
% ZeroX has been modified so as to avoid the following error "Error using griddedInterpolant
% The grid vectors are not strictly monotonic increasing."
%
% The following modifications are based off of this MatLab forum answer:
% https://www.mathworks.com/matlabcentral/answers/283718-grid-vectors-not-strictly-monotonic-increasing
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns Approximate Zero-Crossing Indices Of Argument Vector
zxidx = zci(y);
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
end