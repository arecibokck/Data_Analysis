%Only Use High Resolution data

%read CSV file
filename = 'E:\Master Thesis\Data\MOTLaser_MonitorPD.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f';
fileID = fopen(filename,'r');
csvdat = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

Data = cell2mat(csvdat);
LPower = Data(:,1);
MPDReadings = Data(:,2);

plot(MPDReadings,LPower);
zoom on; grid on;
fTitle  = title(['MOT Laser Monitor PD Response']);
fXLabel = xlabel('PD Readings [mV]'); 
fYLabel = ylabel('Laser Power[mW]');
set( gca                       , ...
    'FontName'   , 'Calibri Light' );
set([fTitle, fXLabel, fYLabel], ...
    'FontName'   , 'Calibri Light');
set([fXLabel, fYLabel]  , ...
    'FontSize'   , 16          );
set( fTitle                    , ...
    'FontSize'   , 16          );%, ...
    %'FontWeight' , 'bold' 
print(gcf,[filename(1:find(filename == '\' , 1, 'last')) 'MOTMPDReadings.png'],'-dpng','-r300'); 