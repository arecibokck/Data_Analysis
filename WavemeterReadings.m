%Only Use High Resolution data

%read CSV file
filename = 'E:\Master Thesis\Data\Optical Path Length Measurements\LaserFrequencyScanforPLD\OPLL active for one arm\HDT3 latched\scanup.csv';
delimiter = ',';
startRow = 2;
% formatSpec = '%f';
formatSpec = '%f%f%f';
fileID = fopen(filename,'r');
WMTrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

Data = cell2mat(WMTrace);
Scan = Data(:,2);
Times = (Data(:,1) - Data(1,1)).*1e-03;
% n1 = length(Scan) - 2.445e04; %HDT1L
% n2 = length(Scan) - 1000; %HDT1L
% n1 = 1; %HDT1R
% n2 = length(Scan) - 115; %HDT1R
% n1 = 1.12e04; %HDT3L
% n2 = length(Scan) - 3.09e04; %HDT3L
% n1 = 8.33e03; %HDT3R
% n2 = length(Scan)-55; %HDT3R
% n1 = 500; %HDT1LR&3LR
% n2 = length(Scan)-780; %HDT1LR&3LR
n1 = 330;
n2 = 352;
Scan = Scan(n1:n2);
Times = Times(n1:n2) - Times(n1);
Scan = (299792458)./ (Scan*1e-09);
%% Plot Readings
figure(1)
hold on
zoom on
grid on
plot(Times, Scan);
range = max(Scan) - min(Scan);
n = find(Scan==max(Scan));
Scan = Scan(n:end);
Times = Times(n:end);
rate = mean(abs(diff(Scan)) ./ diff(Times));
hold off
xlabel('time [s]')
ylabel('Laser Frequency [THz]')
%title(['IFL Laser Piezo Scan (-3 V to 3 V @ 1 Hz):' ' \Delta\nu = ' num2str(range * 1e03, '%.f') ' MHz' ])
title(['Matisse CS Reference Cell Scan (@ ' num2str(sprintf('%0.1fe%+03d',rate*10^-6,6)) ' Hz/s):' ' \Delta\nu = ' num2str(sprintf('%0.1fe%+03d',range*10^-9,9)) ' Hz' ])
%'Range = ' num2str(range) ' Ghz'
xlim([0 max(Times)])
ytix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(ytix/1e12))
% xtix=get(gca,'xtick')';
% set(gca,'xticklabel',num2str(xtix))
print(gcf,[filename(1:find(filename == '\' , 1, 'last')) 'scan.png'],'-dpng','-r300');   