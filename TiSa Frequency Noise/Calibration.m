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
truncated_V_PD = smooth(normV_PD(ll:up));
truncated_time = time(ll:up);
truncated_scan = smooth(scan(ll:up));

%Detect peaks (tops of fringes), extracting both voltage and corresponding time
pks = findpeaksx(truncated_time,truncated_V_PD, 0.001, 0.2, 2, 2, 3);
top_of_fringe_Vs = pks(:,3);
top_of_fringe_times = pks(:,2);
top_of_fringe_V = max(top_of_fringe_Vs);
top_of_fringe_time = top_of_fringe_times(find(top_of_fringe_Vs == max(top_of_fringe_Vs)));

%Numerical values of required parameters
fringe_amplitude = (max(V_PD(:)) - min(V_PD(:)))*1E3;

[minValue_t,closestIndex_t] = min(abs(truncated_time - (top_of_fringe_time - (0.4 * mean(diff(top_of_fringe_times))))));
Index_T = find(truncated_time == top_of_fringe_time);
t_interval = closestIndex_t:Index_T;
rising_edge = truncated_V_PD(t_interval);
[minValue_sof_10,closestIndex_sof_10] = min(abs(rising_edge - (0.1 * top_of_fringe_V)));
[minValue_sof_80,closestIndex_sof_80] = min(abs(rising_edge - (0.8 * top_of_fringe_V)));
[minValue_sof_90,closestIndex_sof_90] = min(abs(rising_edge - (0.9 * top_of_fringe_V)));
rise_time = (truncated_time(t_interval(closestIndex_sof_90)) - truncated_time(t_interval(closestIndex_sof_10)))*1E6;

FSR = 1.5;
t2t1_interval = mean(diff(top_of_fringe_times)*1E6);
conversion_factor = 0.8*((fringe_amplitude*1E-3)/(rise_time*1E-6))*((t2t1_interval*1E-6)/(FSR*1E9));

%Plot
hold on;
plot(truncated_time, truncated_scan);
plot(truncated_time, truncated_V_PD);
plot(truncated_time(t_interval(closestIndex_sof_80)), rising_edge(closestIndex_sof_80), 'bo');
text(truncated_time(t_interval(closestIndex_sof_80)) - 8E-5, rising_edge(closestIndex_sof_80), 'Lock Point -------->');
line([truncated_time(1), truncated_time(end)], [rising_edge(closestIndex_sof_10), rising_edge(closestIndex_sof_10)], 'LineStyle', '--');
line([truncated_time(1), truncated_time(end)], [rising_edge(closestIndex_sof_90), rising_edge(closestIndex_sof_90)], 'LineStyle', '--');
line([truncated_time(t_interval(closestIndex_sof_10)), truncated_time(t_interval(closestIndex_sof_10))],...
    [0, truncated_scan(end)], 'LineStyle', '--');
line([truncated_time(t_interval(closestIndex_sof_90)), truncated_time(t_interval(closestIndex_sof_90))],...
    [0, truncated_scan(end)], 'LineStyle', '--');
text(min(truncated_time)+3E-05, max(truncated_scan)-0.6, ['U_{peak} = ' num2str(fringe_amplitude) ' mV']);
text(min(truncated_time)+3E-05, max(truncated_scan)-0.8, ['T_{rise} = ' num2str(rise_time) ' us']);
text(min(truncated_time)+3E-05, max(truncated_scan)-1.0, ['FSR = ' num2str(FSR) ' Ghz']);
text(min(truncated_time)+3E-05, max(truncated_scan)-1.2, ['t2-t1 = ' num2str(t2t1_interval) ' us' ]);
text(min(truncated_time)+3E-05, max(truncated_scan)-1.4, ['\DeltaU / \Delta \nu = ' num2str(conversion_factor) ' V/Hz' ]);
hold off;
axis tight;
title('Side-Of-Fringe Lock for Fabry-Perot Etalon');
legend('Etalon Piezo Scan','Transmission through Etalon (Normalized)', 'Location', 'northwest');
xlabel('Time(s)') 
ylabel('Voltage(V)')