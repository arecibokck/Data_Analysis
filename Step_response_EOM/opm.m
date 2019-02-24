filename = 'C:\Users\KARTHIK KC\Documents\MATLAB\transmission.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
OPMtrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

voltage = OPMtrace{:, 1};
vals = OPMtrace{:, 2};
power = (vals - min(vals)) / (max(vals) - min(vals));

plot(voltage, power);

axis tight;
title('Transmission through EOM with crossed polarizers');
legend('Power', 'Location', 'northwest');
xlabel('U_{EOM}(V)') 
ylabel('Normalized Power Reading')
clearvars