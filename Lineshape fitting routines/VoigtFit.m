%Only Use High Resolution data

%read CSV file
filename = 'C:\Users\chakar\Documents\MATLAB\SDT2D\adwin\trunk\Scripts\+Karthik\Data Analysis\Lineshape fitting routines\Line_MOT.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f';
fileID = fopen(filename,'r');
LineData = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

freq = LineData{1};
pd = LineData{2};

loglog(freq, pd);
zoom on; grid on;
fTitle  = title(['']);
fXLabel = xlabel('Freq [Hz]'); 
fYLabel = ylabel('Power Density [dB/Hz]');
xlim([100 1E+06])