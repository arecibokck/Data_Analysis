%Only Use High Resolution data

%read CSV file
% filename = 'C:\Users\chakar\Desktop\2020-03-03_DQSIM_TiSa_Linewidth Measurements\PSD_latticelaser.csv';
filename = 'C:\Users\chakar\Desktop\2020-03-03_DQSIM_TiSa_Linewidth Measurements\PSD_MOTlaser.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f';
fileID = fopen(filename,'r');
PSDData = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
%%

freq = PSDData{1};
psd = PSDData{2};

loglog(freq, psd);
zoom on; grid on;
fTitle  = title(['Freq Noise']);
fXLabel = xlabel('Freq [Hz]'); 
fYLabel = ylabel('Noise Spectral Density [Hz/sqrt(Hz)]');
xlim([100 1E+06])
set( gca                       , ...
    'FontName'   , 'Calibri Light' );
set([fTitle, fXLabel, fYLabel], ...
    'FontName'   , 'Calibri Light');
set([fXLabel, fYLabel]  , ...
    'FontSize'   , 16          );
set( fTitle                    , ...
    'FontSize'   , 16          );%, ...
    %'FontWeight' , 'bold' 