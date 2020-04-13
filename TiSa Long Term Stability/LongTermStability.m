%read CSV file
filename = 'E:\Master Thesis\Data\TiSa Laser Long-term Frequency Stability Measurement\20190922\Matisse CR\scan.csv';
delimiter = ',';
startRow = 2;
% formatSpec = '%f';
formatSpec = '%f%f%f';
fileID = fopen(filename,'r');
WMTrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

Data = cell2mat(WMTrace);
Time = Data(:,1);
w = Data(:,2);
%w = w(~isnan(w)); 
Temp = Data(:,3);
f = (299792458)./ (w*1e-09);

% f = f(1:1e03:end); %downsample(f,1e03);
% f = smoothdata(f, 'gaussian', 30);
t = (linspace(0,max(Time)*1e-03,length(f)))./3600; 
% t = (linspace(0,max(Time)*1e-03,length(f)))./60;
totaltime = (Time(end)*1e-03)/3600;

%% Allan Deviation
tau0 = 5;
avar = [];
tau = [];
x=[];
n=length(f);
jj=floor(log((n-1)/3)/log(2));
for j=0:jj
    m=2^j;
    tau(j+1)=m*tau0;
    D=zeros(1,n-m+1);
    for i=1:n-m+1
        D(i)=sum(f(i:i+m-1))/m;
    end
    avar(j+1)=sqrt(0.5*mean((diff(D(1:m:n-m+1)).^2)));
end
%% Plot Frequency & Temperature vs Time
figure(1)
% f = smoothdata(f, 'gaussian', 30);
% yyaxis left
% set(gca, 'YColor', [0, 0, 0]);
% plot(t, f, 'LineWidth',1.5);
plot(t, (f-f(1))*1E-6, 'LineWidth',1.5);
% plot(t, (f-mean(f))*1E-6, 'LineWidth',1.5);
grid on
xlim([0 max(t)])
% ytix=get(gca,'ytick')';
% set(gca,'yticklabel',num2str(ytix/1e12, '%0.6f'))
% fTitle  = title(['Matisse CR Long Term Stability: Measured for ' num2str(totaltime, '%0.f') ' Hrs']);
fTitle  = title(['Matisse CR Long Term Stability']);
% fTitle  = title(['Matisse CS locked to Transfer Cavity: Short-term Stability']);
% fTitle  = title(['Matisse CS locked to Transfer Cavity: Long-term Stability']);
fXLabel = xlabel('time [Hrs]'); 
% fXLabel = xlabel('time [min]'); 
% fXLabel = xlabel('time [s]'); 
fY1Label = ylabel('Drift of Laser Frequency [MHz]');
% yyaxis right
% set(gca, 'YColor', [0, 0, 0]);
% plot(t, Temp, 'LineWidth',1.5, 'Color', [0.4940, 0.1840, 0.5560]); 
% fY2Label = ylabel(['Temperature[' char(176) 'C]']);
fLegend = legend(...
          'Laser Frequency', 'Temperature', 'Location','northeast');
          %['Total Drift = ' num2str(sprintf('%0.1fe%+03d',drift*10^-6,6)) ' Hz'], 'Temperature', 'Location','northeast');
set( gca                       , ...
    'FontName'   , 'Calibri Light' );
set([fTitle, fXLabel, fY1Label], ...
    'FontName'   , 'Calibri Light');
% set([fTitle, fXLabel, fY1Label, fY2Label], ...
%     'FontName'   , 'Calibri Light');
set([fLegend, gca]             , ...
    'FontSize'   , 16           );
set([fXLabel, fY1Label]  , ...
    'FontSize'   , 14          );
% set([fXLabel, fY1Label, fY2Label]  , ...
%     'FontSize'   , 14          );
set( fTitle                    , ...
    'FontSize'   , 14          );%, ...
    %'FontWeight' , 'bold' 
print(gcf,[filename(1:find(filename == '\' , 1, 'last')) 'LTS_Drift.png'],'-dpng','-r300'); 

%% Plot Allan Deviation

figure(2)
loglog(tau, avar*1E-6, '-s', 'LineWidth',1.5);
% - adding vertical lines for minute, hour and day
% h1 = vline(60,'b','1 minute');
% h2 = vline(60*60,'r','1 hour');
% h3 = vline(60*60*24,'g','1 day'); 
zoom on; grid on;
fTitle  = title(['Allan Deviation of Laser Frequency']);
fXLabel = xlabel('\tau [s]'); 
fYLabel = ylabel('\sigma_y(\tau) [MHz]');
set( gca                       , ...
    'FontName'   , 'Calibri Light' );
set([fTitle, fXLabel, fYLabel], ...
    'FontName'   , 'Calibri Light');
set([fXLabel, fYLabel]  , ...
    'FontSize'   , 16          );
set( fTitle                    , ...
    'FontSize'   , 16          );%, ...
    %'FontWeight' , 'bold' 
print(gcf,[filename(1:find(filename == '\' , 1, 'last')) 'MatisseCR_AllanDev.png'],'-dpng','-r300'); 