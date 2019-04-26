%read CSV file
filename = 'scope_1.csv';
delimiter = ',';
startRow = 3;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
OPMtrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

%Extract piezo scan (the rising edge) and the corresponding photodector signal
scan_channel = 3;
V_PD_signal_channel = 2;
time = OPMtrace{1};
time = time - min(time);
scan = OPMtrace{scan_channel} - min(OPMtrace{scan_channel});
V_PD = OPMtrace{V_PD_signal_channel};
normV_PD = (V_PD - min(V_PD(:))) ./ (max(V_PD(:)) - min(V_PD(:)));
ll = find(scan == min(scan));
up = find(scan == max(scan));
truncated_V_PD = normV_PD(ll:up);
truncated_time = time(ll:up);
truncated_scan = scan(ll:up);

%Detect peaks (tops of fringes), extracting both voltage and corresponding time
pks = findpeaksx(truncated_time,truncated_V_PD, 0.001, 0.2, 2, 2, 3);
top_of_fringe_Vs = pks(:,3);
top_of_fringe_times = pks(:,2);
top_of_fringe_V = max(top_of_fringe_Vs);
top_of_fringe_time = top_of_fringe_times(find(top_of_fringe_Vs == max(top_of_fringe_Vs)));

%Numerical values of required parameters
fringe_amplitude = (max(V_PD(:)) - min(V_PD(:)))*1E3;
rise_time = 14.39;
FSR = 1.5;
t2t1_interval = mean(diff(top_of_fringe_times)*1E6);

%Plot
hold on;
plot(truncated_time, truncated_scan);
plot(truncated_time, truncated_V_PD);
text(min(truncated_time)+3E-05, max(truncated_scan)-0.6, ['U_{peak} = ' num2str(fringe_amplitude) ' mV']);
text(min(truncated_time)+3E-05, max(truncated_scan)-0.8, ['T_{rise} = ' num2str(rise_time) ' us']);
text(min(truncated_time)+3E-05, max(truncated_scan)-1.0, ['FSR = ' num2str(FSR) ' Ghz']);
text(min(truncated_time)+3E-05, max(truncated_scan)-1.2, ['t2-t1 = ' num2str(t2t1_interval) ' us' ]);
hold off;
axis tight;
title('Side-Of-Fringe Lock for Fabry-Perot Etalon');
legend('Etalon Piezo Scan','Transmission through Etalon (Normalized)', 'Location', 'northwest');
xlabel('Time(s)') 
ylabel('Voltage(V)')