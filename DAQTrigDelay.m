%Only Use High Resolution data

%read CSV file
filename = 'E:\Master Thesis\Code\Data Analysis\DAQ Trigger Delay\scope_comb.csv';
delimiter = ',';
startRow = 3;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
ScopeTrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

%Plot
figure(1);
hold on;
time = ScopeTrace{1} - ScopeTrace{1}(1);
temp = ScopeTrace{2};
temp(temp <= 1) = 0;
temp(temp >  1) = 1;
in = find(temp, 1, 'first');
trigTime = time(in);
plot(time, ScopeTrace{2});
Delay = [];
for i = 3:length(ScopeTrace)-1
    temp = ScopeTrace{i};
    temp(temp <= 1) = 0;
    temp(temp >  1) = 1;
    in = find(temp, 1, 'first');
    Delay = [Delay (time(in)-trigTime)];
    plot(time, ScopeTrace{i});
end
xlabel('time [s]')
ylabel('Trigger State')
% export_fig('C:\Users\chakar\Desktop\TrigState', '-png');

hold off

figure(2)
hold on
plot(Delay)
mn = mean(Delay);
er = std(Delay);
line([1, length(Delay)], [mn, mn], 'LineStyle', '--', 'Color',[.5 .4 .7])
text(5, mn+5e-04, ['Mean = ' num2str(mn*1e+03) ' ± ' num2str(sprintf('%.3f',er*1e+03)) ' ms']);
ytix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(ytix*1e+03))
xlabel('Run')
ylabel('Delay at Start (ms)')
% export_fig('C:\Users\chakar\Desktop\StartDelay', '-png');

temp = ScopeTrace{2};
temp(temp <= 1) = 0;
temp(temp >  1) = 1;
in = find(temp, 1, 'last');
trigTime = time(in);
Delay = [];
for i = 3:length(ScopeTrace)-1
    temp = ScopeTrace{i};
    temp(temp <= 1) = 0;
    temp(temp >  1) = 1;
    in = find(temp, 1, 'last');
    Delay = [Delay (time(in)-trigTime)];
end

figure(3)
hold on
plot(Delay)
mn = mean(Delay);
er = std(Delay);
line([1, length(Delay)], [mn, mn], 'LineStyle', '--', 'Color',[.5 .4 .7])
text(5, mn+5e-04, ['Mean = ' num2str(mn*1e+03) ' ± ' num2str(sprintf('%.3f',er*1e+03)) ' ms']);
ytix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(ytix*1e+03))
xlabel('Run')
ylabel('Delay at Stop (ms)')
% export_fig('C:\Users\chakar\Desktop\StopDelay', '-png');
