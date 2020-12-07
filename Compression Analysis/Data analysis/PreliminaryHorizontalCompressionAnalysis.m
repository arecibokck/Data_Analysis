
% load('2020-02-20T215858_seq_atomDynamics.mat');
% meas = meas20200220T215858;

load('C:\Users\chakar\Documents\MATLAB\SDT2D\adwin\trunk\Scripts\+Karthik\Data Analysis\Compression Analysis\20200410\2020-04-10T045138_seq_scanHorizontalCompression.mat')
meas = meas20200410T045138;
% method = 'CrossCorrelation';
method = 'WeightedAveraging';

mdiffX = [];
errX = [];
mdiffY = [];
errY = [];

% figure('Position', [70, 500, 1800, 400])
figure('Position', [70, 70, 1200, 900])
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
% m = 21;
% n = 6;
% m = size(meas.runData.usedScanParameters.probeTime, 1);
m = size(meas.runData.usedScanParameters.WaitTimeHorizontal,1);
% m = size(meas.runData.usedScanParameters.VDTPowerDuringHorizontalCompression,1);
n = size(meas.runData.usedScanParameters.WaitTimeHorizontal,2);
% n = size(meas.runData.usedScanParameters.VDTPowerDuringHorizontalCompression,2);
% pbtimes = meas.runData.usedScanParameters.probeTime(1:m);
pbtimes = meas.runData.usedScanParameters.WaitTimeHorizontal(1:m);
% VDTpowers = meas.runData.usedScanParameters.VDTPowerDuringHorizontalCompression(1:m);
Images = cell (m, n);
for i = 1:m
    c = 1;
    for j = 0:m:size(meas.runData.images,1)-1
        Images{i,c} = squeeze(meas.runData.images(i+j,:,:,:));
        c = c + 1;
    end
end

%% Tracking the center of mass of atom cloud
for i = 1:m
    diffinpixelsX =[];
    diffinpixelsY = [];
    for j = 1 : n
        img = Images{i,j};
        
        backgroundImage = squeeze(img(3,:,:));
        
        firstImage = squeeze(img(1,:,:)) - backgroundImage;
        firstImage = firstImage - min (firstImage (:));
        firstImage = firstImage / max (firstImage (:));

        secondImage = squeeze(img(2,:,:)) - backgroundImage;
        secondImage = secondImage - min (secondImage (:));
        secondImage = secondImage / max (secondImage (:));
        
        switch method
            
            case 'CrossCorrelation'
                usfac = 100;
                [output, Greg] = dftregistration(fft2(firstImage),fft2(secondImage),usfac);
                diffinpixelsX(end+1) = abs(output(4));
                diffinpixelsY(end+1) = abs(output(3));
        
            case 'WeightedAveraging'
                fi = Wavg(firstImage);
                si = Wavg(secondImage);
                diffinpixelsX(end+1) = abs(fi(1) - si(1));
                diffinpixelsY(end+1) = abs(fi(2) - si(2));   
        end
        
        subplot(2,2,1), imagesc(firstImage)
        title('Pre-Compression')
        
        subplot(2,2,2), imagesc(secondImage)
        title('Post-Compression')
        
    end
    mx = mean(diffinpixelsX);
    my = mean(diffinpixelsY);
%   mdiffX(end+1) = (1/sqrt(2)*(mx - my))/1.69;
%   errX(end+1) = std((1/sqrt(2)*(diffinpixelsX - diffinpixelsY)))/sqrt(size(diffinpixelsX,2));
    mdiffX(end+1) = mx/1.69;
    errX(end+1) = std(diffinpixelsX/1.69)/sqrt(size(diffinpixelsX,2));
%   mdiffY(end+1) = (1/sqrt(2)*(mx + my))/1.69;
%   errY(end+1) = std((1/sqrt(2)*(diffinpixelsX + diffinpixelsY)))/sqrt(size(diffinpixelsY,2));
    mdiffY(end+1) = my/1.69;
    errY(end+1) = std(diffinpixelsY/1.69)/sqrt(size(diffinpixelsY,2));
    
end
subplot(2,2,3)
hold on
xshift = plot(pbtimes, mdiffX, 'o--', 'MarkerFaceColor', colours{1}, 'MarkerEdgeColor', colours{1}, 'MarkerSize', 4);
% xshift = plot(VDTpowers, mdiffX, 'o--', 'MarkerFaceColor', colours{1}, 'MarkerEdgeColor', colours{1}, 'MarkerSize', 4);
errorbar (pbtimes, mdiffX, errX, 'Color', colours{1}, 'LineStyle', 'None')
% errorbar (VDTpowers, mdiffX, errX, 'Color', colours{1}, 'LineStyle', 'None')
yshift = plot(pbtimes, mdiffY, 'o--', 'MarkerFaceColor', colours{2}, 'MarkerEdgeColor', colours{2}, 'MarkerSize', 4);
% yshift = plot(VDTpowers, mdiffY, 'o--', 'MarkerFaceColor', colours{2}, 'MarkerEdgeColor', colours{2}, 'MarkerSize', 4);
errorbar (pbtimes, mdiffY, errY, 'Color', colours{2}, 'LineStyle', 'None')
% errorbar (VDTpowers, mdiffY, errY, 'Color', colours{2}, 'LineStyle', 'None')
xlim([0 3001]);
% xlim([0 2]);
xlabel('Probe Time (microseconds)');
% xlabel('VDT power (Watts)');
ylabel('Measured shift (# of Lattice Sites)');
legend ([xshift, yshift], { 'Along image X-axis' , 'Along image Y-axis' }, 'Location', 'northwest');
title('Shift of CoM of atom cloud');

%% Tracking the RMS spread of atom cloud

f_X_spread = [];
err_f_X_spread = [];
f_Y_spread = [];
err_f_Y_spread = [];
g_X_spread = [];
err_g_X_spread = [];
g_Y_spread = [];
err_g_Y_spread = [];

for i = 1:m
    f_xspread = [];
    f_yspread = [];
    g_xspread = [];
    g_yspread = [];
    
    for j = 1 : n
        img = Images{i,j};
        
        backgroundImage = squeeze(img(3,:,:));
        
        firstImage = squeeze(img(1,:,:)) - backgroundImage;
        firstImage = firstImage - min (firstImage (:));
        firstImage = firstImage / max (firstImage (:));
        
        f_xspread(end+1) = std(sum(firstImage, 1)); % std(sum(spdiags(f)));
        f_yspread(end+1) = std(sum(firstImage, 2)); %std(sum(spdiags(rot90(f))));
        
        secondImage = squeeze(img(2,:,:)) - backgroundImage;
        secondImage = secondImage - min (secondImage (:));
        secondImage = secondImage / max (secondImage (:));

        g_xspread(end+1) = std(sum(secondImage, 1)); % std(sum(spdiags(g)));
        g_yspread(end+1) = std(sum(secondImage, 2)); % std(sum(spdiags(rot90(g))));
         
        
    end
    cf = 0.612/1.69;
    f_X_spread(end+1) = mean(f_xspread) * cf;...(mean(f_xspread) - mean(g_xspread))/1.69;
    err_f_X_spread(end+1) = std(f_xspread .* cf)/sqrt(size(f_xspread,2));
    f_Y_spread(end+1) = mean(f_yspread) * cf; ...(mean(f_yspread) - mean(g_yspread))/1.69;
    err_f_Y_spread(end+1) = std(f_yspread .* cf)/sqrt(size(f_yspread,2));
    g_X_spread(end+1) = mean(g_xspread) * cf;...(mean(f_xspread) - mean(g_xspread))/1.69;
    err_g_X_spread(end+1) = std(g_xspread .* cf)/sqrt(size(g_xspread,2));
    g_Y_spread(end+1) = mean(g_yspread) * cf; ...(mean(f_yspread) - mean(g_yspread))/1.69;
    err_g_Y_spread(end+1) = std(g_yspread .* cf)/sqrt(size(g_yspread,2));
end

subplot(2,2,4)
hold on
xsp = plot(pbtimes, (g_X_spread./f_X_spread), 'o--', 'MarkerFaceColor', colours{1}, 'MarkerEdgeColor', colours{1}, 'MarkerSize', 4);
% xsp = plot(VDTpowers, (g_X_spread./f_X_spread), 'o--', 'MarkerFaceColor', colours{1}, 'MarkerEdgeColor', colours{1}, 'MarkerSize', 4);
errorX = sqrt((err_f_X_spread./f_X_spread).^2 + (err_g_X_spread./g_X_spread).^2);
errorbar (pbtimes, (g_X_spread./f_X_spread), errorX, 'Color', colours{1}, 'LineStyle', 'None')
% errorbar (VDTpowers, (g_X_spread./f_X_spread), errorX, 'Color', colours{1}, 'LineStyle', 'None')
ysp = plot(pbtimes, (g_Y_spread./f_Y_spread), 'o--', 'MarkerFaceColor', colours{2}, 'MarkerEdgeColor', colours{2}, 'MarkerSize', 4);
% ysp = plot(VDTpowers, (g_Y_spread./f_Y_spread), 'o--', 'MarkerFaceColor', colours{2}, 'MarkerEdgeColor', colours{2}, 'MarkerSize', 4);
errorY = sqrt((err_f_Y_spread./f_Y_spread).^2 + (err_g_Y_spread./g_Y_spread).^2);
errorbar (pbtimes, (g_Y_spread./f_Y_spread), errorY, 'Color', colours{2}, 'LineStyle', 'None')
% errorbar (VDTpowers, (g_Y_spread./f_Y_spread), errorY, 'Color', colours{2}, 'LineStyle', 'None')
xlim([0 3001]);
xlabel('Probe Time (microseconds)');
% xlabel('VDT Power (Watts)');
ylabel('Spread ratio');
lg = legend ([xsp, ysp], { 'Along image X-axis' , 'Along image Y-axis' }, 'Location', 'northeast');
title(lg, 'Post:Pre');
title('RMS spread of atom cloud');

function output = Wavg(A)
    %CENTEROFMASS Summary of this function goes here
    %   Detailed explanation goes here

    A = squeeze(A); % - remove singleton dimension
    % Method 1, using mean
    x = 1 : size(A, 2); % Columns.
    y = 1 : size(A, 1); % Rows.
    [X, Y] = meshgrid(x, y);
    meanA = mean(A(:));
    centerOfMassX = mean(A(:) .* X(:)) / meanA;
    centerOfMassY = mean(A(:) .* Y(:)) / meanA;
    output = [centerOfMassX, centerOfMassY, meanA];
end