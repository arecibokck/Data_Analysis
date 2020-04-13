%Only Use High Resolution data

%read CSV file
filename = 'E:\Master Thesis\Data\Optical Path Length Measurements\At Experiment\HDT1L\Initial Path Length\beathdt1l.csv';
delimiter = ',';
startRow = 3;
formatSpec = '%f%f';
fileID = fopen(filename,'r');
DMMTrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

Data = cell2mat(DMMTrace);
Time = Data(:,1) - Data(1,1);
Beat = Data(:,2);
% i =  find(Time==15)+40; %HDT1L                                             1; %HDT1R
% j =  find(Time==95); %HDT1L                                   find(Time==70); %HDT1R
% t = Time(i:j) - Time(i);
ti = 55; %HDT1L & %HDT1R
tf = 65; %HDT1L & %HDT1R 
% t = Time;
% Beat = Beat(~isnan(Beat)); %HDT3R
% t = Time(1:length(Beat)); %HDT3R
% ti = 65; %40; %HDT3LR
% tf = 85; % 90;full scan for HDT3LR
% ti = 60;   %35;
% tf = 85;   %85;full scan for HDT1LR
t = Time(find(Time==ti):find(Time==tf)) - Time(find(Time==ti));
Beat = detrend(Beat(find(Time==ti):find(Time==tf)) - Beat(find(Time==ti)));
%% Plot Readings
figure(1)
hold on
% i = find(Scan == min(Scan), 1, 'first');
% Beat = Beat([i:i+500000]);
% Times = Times([i:i+500000]) - Times(i);
% Scan = Scan([i:i+500000]);
axis tight
plot(t, Beat)
% plot(Times, Scan)
hold off
fTitle  = title ('Beat Signal');
fXLabel = xlabel('time [ms]'); 
fYLabel = ylabel('Voltage [V]');
% xlim([0 Times(end:end)])
set( gca                       , ...
    'FontName'   , 'Calibri Light' );
set([fTitle, fXLabel, fYLabel], ...
    'FontName'   , 'Calibri Light');
set([fXLabel, fYLabel]  , ...
    'FontSize'   , 14          );
set( fTitle                    , ...
    'FontSize'   , 16          );%, ...
    %'FontWeight' , 'bold'      
% grid
%print(gcf,'5m_fiber.png','-dpng','-r300');
%% Fitting a Sinusoid
% y = detrend(downsample(Beat,1e01));
% y = downsample(Beat,1e01);
% y = smooth(Beat);
y = Beat;
x = linspace(min(t),max(t),length(y)); 
% y = y(300:end-100);
% x = linspace(0.2990*1e03,1.899*1e03,length(y)); 
yu = max(y);
yl = min(y);
yr = (yu-yl);                                                    % Range of ‘y’
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[1 0]) <= 0);                          % Find zero-crossings
per = 2*mean(diff(zx));                                          % Estimate period
ym = mean(y);                                                    % Estimate offset
fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);      % Function to fit
fcn = @(b) sum((fit(b,x)' - y).^2);                              % Least-Squares cost function
s = fminsearch(fcn, [yr;  per;  -1;  ym]);                       % Minimise Least-Squares
amplitude = s(1);  %(in units of y)
period = s(2);  %(in units of x)
phase = s(3);  %(phase is s(2)/(2*s(3)) in units of x)
offset = s(4);  %(in units of y)
xp = linspace(min(x),max(x), 1000);
figure(2)
dataline = line(t,Beat); 
fitline = line(xp, fit(s,xp));
set(dataline                      , ...
  'LineStyle'       , 'none'      , ...  
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 5           , ...
  'MarkerEdgeColor' , 'none'      , ...
  'MarkerFaceColor' , [0.75 0.75 1] );
set(fitline                       , ...
  'LineStyle'       , '--'        , ...
  'LineWidth'       , 1.5        , ...
  'Color'           , [0.8500, 0.3250, 0.0980] );
fTitle  = title ('Beat Signal');
fXLabel = xlabel('time [s]'); 
fYLabel = ylabel('Voltage [V]');
fLegend = legend( ...
  [fitline, dataline], ...
  ['Fit, \nu_b = ' num2str(1/period)] , ...
  'Data');
set( gca                       , ...
    'FontName'   , 'Calibri Light' );
set([fTitle, fXLabel, fYLabel], ...
    'FontName'   , 'Calibri Light');
set([fLegend, gca]             , ...
    'FontSize'   , 10           );
set([fXLabel, fYLabel]  , ...
    'FontSize'   , 14          );
set( fTitle                    , ...
    'FontSize'   , 16          );%, ...
    %'FontWeight' , 'bold'      
grid
% filename(1:find(filename == '\' , 1, 'last'));
% dataset = split(filename,"\");
% dat = split(dataset(end), '.');
% print(gcf, [filename(1:find(filename == '\' , 1, 'last')) dat{1} '.png'],'-dpng','-r300');   
%% Calculate and Plot FFT
% Beat = smooth(Beat);
% Times = linspace(0,max(t)*1e03,length(Beat)); 
% % f = 0.1;   % frequency Hz
% % L = length(DMMTrace{1});
% % Times = 0:0.5:L; % s
% % X = sin(2 * pi * f * Times)';
% 
% X = mean(Beat,2);
% X = X-mean(X);
% WindowedX = X.*hann(length(X));
% 
% figure(3) 
% subplot(1,3,1)
% plot(Times, X, 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth',1.5)
% xlabel('Time [ms]')
% ylabel('Voltage [V]')
% title('Average Voltage')
% xlim([0 Times(end:end)])
% grid on
% 
% 
% subplot(1,3,2)
% hold off
% plot(Times,WindowedX, 'Color', [0, 0.4470, 0.7410]	, 'LineWidth',1.5)
% xlabel('time [ms]')
% ylabel('Voltage [V]')
% title('Apodization (Window) Applied')
% xlim([0 Times(end:end)])
% hold on
% grid on
% subplot(1,3,3)
% 
% T = 1e-03;              % Sampling period       
% Fs = 1/T;               % Sampling frequency             
% L = length(WindowedX);  % Length of signal
% t = (0:L-1)*T;          % Time vector
% 
% Y = fft(WindowedX);
% 
% %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% %Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes will not be exactly where expected, perhaps because of the noise. On average, longer signals produce better frequency approximations.
% 
% f = Fs * (0:(L/2))/L;
% %f = Fs/2*(0:(L/2));
% hold off
% %plot(f,P1,'.-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth',1.5)
% semilogy(f,P1)
% axis([0  10    ylim])
% hold on
% title('Single-Sided Amplitude Spectrum')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% xlim([0 20])
% grid on