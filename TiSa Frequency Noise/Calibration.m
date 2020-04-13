%Only Use High Resolution data

%read CSV file
filename = 'E:\Master Thesis\Data\Fabry-Perot Resonator\longer.csv';
delimiter = ',';
startRow = 3;
%formatSpec = '%f%f%f%[^\n\r]';
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
OPMtrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
%%
%Extract piezo scan (the rising edge) and the corresponding photodector signal
V_PD_signal_channel = 2;
% scan_channel = 4;
time = OPMtrace{1};
time = time - min(time);
% scan = smooth(OPMtrace{scan_channel} - min(OPMtrace{scan_channel}));
V_PD = OPMtrace{V_PD_signal_channel};
norm_V_PD = (V_PD - min(V_PD(:))) ./ (max(V_PD(:)) - min(V_PD(:)));
% ll = find(scan == min(scan));
% up = find(scan == max(scan));
% ll = find(time == 0.2);
% up = find(time == 0.3);
% truncated_V_PD = 10*V_PD(ll:up);
% truncated_time = time(ll:up);
% truncated_scan = scan(ll:up);

% [yprime2 params2 resnorm2 residual2] = lorentzfit(time, V_PD,[],[],'3');
% figure; plot(x,y,'b.','LineWidth',2)
% hold on; plot(x,yprime2,'r-','LineWidth',2)

%Detect peaks (tops of fringes), extracting both voltage and corresponding time
% pks = findpeaksx(truncated_time,truncated_V_PD, 0.001, 0.2, 2, 2, 3);
% pks = findpeaksx(time,norm_V_PD, 0.001, 0.2, 2, 2, 3);
% top_of_fringe_Vs = pks(:,3);
% top_of_fringe_times = pks(:,2);
% top_of_fringe_V = max(top_of_fringe_Vs);
% top_of_fringe_time = top_of_fringe_times(find(top_of_fringe_Vs == max(top_of_fringe_Vs)));
% 
% % Numerical values of required parameters
% fringe_amplitude = (max(V_PD(:)) - min(V_PD(:)))*1E3;
% 
% [minValue_t,closestIndex_t] = min(abs(truncated_time - (top_of_fringe_time - (0.4 * mean(diff(top_of_fringe_times))))));
% Index_T = find(truncated_time == top_of_fringe_time);
% t_interval = closestIndex_t:Index_T;
% rising_edge = truncated_V_PD(t_interval);
% [minValue_sof_10,closestIndex_sof_10] = min(abs(rising_edge - (0.1 * top_of_fringe_V)));
% [minValue_sof_80,closestIndex_sof_80] = min(abs(rising_edge - (0.286 * top_of_fringe_V)));
% [minValue_sof_90,closestIndex_sof_90] = min(abs(rising_edge - (0.9 * top_of_fringe_V)));
% rise_time = (truncated_time(t_interval(closestIndex_sof_90)) - truncated_time(t_interval(closestIndex_sof_10)))*1E6;
% t2t1_interval = mean(diff(top_of_fringe_times)*1E6);
% 
% FSR = 1.5;
% conversion_factor = 0.8*((fringe_amplitude*1E-3)/(rise_time*1E-6))*((t2t1_interval*1E-6)/(FSR*1E9));

%Plot
figure(1);
hold on;
% plot(time, scan);
plot(time, V_PD);
% plot(truncated_time, truncated_scan);
% plot(truncated_time, truncated_V_PD);
% V_PD = detrend(OPMtrace{3});
% truncated_V_PD = V_PD(ll:up);
% plot(truncated_time, truncated_V_PD);
% plot(truncated_time(t_interval(closestIndex_sof_80)), rising_edge(closestIndex_sof_80)+0.2386, 'bo');
% text(truncated_time(t_interval(closestIndex_sof_80)) + 0.5E-4, rising_edge(closestIndex_sof_80)+0.2386, ' <-------- Lock Point (0.7 V)');
% line([truncated_time(1), truncated_time(end)], [rising_edge(closestIndex_sof_10), rising_edge(closestIndex_sof_10)], 'LineStyle', '--');
% line([truncated_time(1), truncated_time(end)], [rising_edge(closestIndex_sof_90), rising_edge(closestIndex_sof_90)], 'LineStyle', '--');
% line([truncated_time(t_interval(closestIndex_sof_10)), truncated_time(t_interval(closestIndex_sof_10))],...
%     [0, truncated_scan(end)], 'LineStyle', '--');
% line([truncated_time(t_interval(closestIndex_sof_90)), truncated_time(t_interval(closestIndex_sof_90))],...
%     [0, truncated_scan(end)], 'LineStyle', '--');
% text(min(truncated_time)+3E-03, max(truncated_scan)-1.6, ['U_{peak} = ' num2str(fringe_amplitude) ' mV']);
% text(min(truncated_time)+3E-03, max(truncated_scan)-2.0, ['T_{rise} = ' num2str(rise_time) ' us']);
% text(min(truncated_time)+3E-03, max(truncated_scan)-2.4, ['FSR = ' num2str(FSR) ' Ghz']);
% text(min(truncated_time)+3E-03, max(truncated_scan)-2.8, ['t2-t1 = ' num2str(t2t1_interval) ' us' ]);
% text(min(truncated_time)+3E-03, max(truncated_scan)-3.2, ['\DeltaU / \Delta \nu = ' num2str(conversion_factor) ' V/Hz' ]);
hold off;
axis tight;
title('Side-Of-Fringe Lock for Fabry-Perot Cavity');
legend('Cavity Piezo Scan','Transmission through Cavity', 'Location', 'northwest');
xlabel('Time(s)') 
ylabel('Voltage(V)')
print(1, 'E:\Master Thesis\Data\TiSa Laser Frequency Noise\Matisse CS\190915_SO_LPF\Locked to 866\Scope\ScopeTrace.png', '-dpng','-r300');