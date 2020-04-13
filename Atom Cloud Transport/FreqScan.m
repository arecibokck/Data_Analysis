%Only Use High Resolution data

%read CSV file
format long
filename = 'E:\Master Thesis\Data\Optical Path Length Measurements\2019-12-05_LaserFrequencyScanforPLD\OPLL active for one arm\HDT1 latched\05.12.2019, 22.13,  346,187849 THz.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f';
fileID = fopen(filename,'r');
WMTrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

Data = cell2mat(WMTrace);
Scan = (299792458 ./ Data(:,2)) *1e09;
Times = (Data(:,1) - Data(1,1)).*1e-03;

%% Plot Readings
start = 600;
stop = 670;
t = Times(start:stop) - Times(start) ;
s = Scan(start:stop);

figure(1)
hold on
zoom on
grid on
plot(t,s);
range = max(s) - min(s)
rate = mean(abs(diff(s)) ./ diff(t))
hold off
xlabel('time [s]')
ylabel('Laser Frequency [THz]')
title(['Matisse CS Reference Cell Scan (@ ' num2str(sprintf('%0.1fe%+03d',rate*10^-6,6)) ' Hz/s):' ' \Delta\nu = ' num2str(sprintf('%0.1fe%+03d',range*10^-9,9)) ' Hz' ])