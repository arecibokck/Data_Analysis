% load('\\rua.iap.uni-bonn.de\users\Karthik\nobackup_LabData\AtomCloudCenterOfMass\2019-11-19T192230_seq_transportAtoms.mat')
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-11-19T192230_seq_transportAtoms.mat')
meas = meas20191119T192230;
mdiff = [];
err = [];
figure('Position', [70, 70, 1200, 900])
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
n = 5;
SortedImages = cell (n, 10);
for i = 1:n
    c = 1;
    for j = 0:n:size(meas.runData.images, 1)-1
       SortedImages{i,c} = squeeze(meas.runData.images(i+j,:,:,:));
       c = c + 1;
    end
end

frame = 1;
for i = 1:n
diffinpixelsX =[];
diffinpixelsY = [];
for j = 1 : 10

img = SortedImages{i,j};
img1 = squeeze(img(1,:,:));
img1 = img1 - min (img1 (:));
img1 = img1 / max (img1 (:));

img2 = squeeze(img(2,:,:));
img2 = img2 - min (img2 (:));
img2 = img2 / max (img2 (:));

autocorr = normxcorr2_general(img1, img1);
autocorr(1:10,:) = NaN;
autocorr(:,1:10) = NaN;
autocorr(:,end-10:end) = NaN;
autocorr(end-10:end,:) = NaN;
[acypeak,  acxpeak] = find (autocorr == max (autocorr (:)));

if acxpeak > size(img1,1)
    acxpeak = acxpeak - size(img1,1);
end
if acypeak > size(img1,2)
    acypeak = acypeak - size(img1,2);
end
    

crosscorr = normxcorr2_general(img2, img1);
crosscorr(1:10,:) = NaN;
crosscorr(:,1:10) = NaN;
crosscorr(:,end-10:end) = NaN;
crosscorr(end-10:end,:) = NaN;
[ccypeak,  ccxpeak] = find (crosscorr == max (crosscorr (:)));

if ccxpeak > size(img2,1) 
    ccxpeak = ccxpeak - size(img2,1);
end

if ccypeak > size(img2,2)
    ccypeak = ccypeak - size(img2,2);
end

subplot(2,2,1), imagesc(img1)
% hold on
% plot(0.5*(acxpeak+1),0.5*(acypeak+1), '+', 'Color', [1 1 1], 'MarkerSize', 10, 'LineWidth', 2);
% hold off
colorbar
title('First Image')

subplot(2,2,2), imagesc(img2)
% hold on
% plot(0.5*(ccxpeak+1),0.5*(ccypeak+1), '+', 'Color', [1 1 1], 'MarkerSize', 10, 'LineWidth', 2);
% hold off
colorbar
title('Second Image')

subplot(2,2,3)
ac = surf(autocorr, 'FaceColor' , colours{1} , 'FaceAlpha' , 0.9, 'EdgeColor' , 'none' );
hold on
cc = surf(crosscorr, 'FaceColor' , colours{3} , 'FaceAlpha' , 0.9, 'EdgeColor' , 'none' );
legend ([ac, cc], { 'First Image Auto-Correlation' , 'Cross-Correlation' }, 'Location', 'best');
title('Feature matching through Spatial Intensity Correlation')
% shading flat
hold off

diffinpixelsX(end+1) = abs(ccxpeak - acxpeak);
diffinpixelsY(end+1) = abs(ccypeak - acypeak);
subplot(2,2,4)
hold on
xlim([1 5]);
xlabel('Index');
ylabel('Along direction of shift (# of Lattice Sites)');
title('Shift of atom cloud')
F (frame) = getframe (gcf);
drawnow
frame = frame+1;
end
mdiff(end+1) = mean(diffinpixelsY)/1.5;
err(end+1) = std(diffinpixelsY);
subplot(2,2,4), plot(mdiff, 'o', 'MarkerFaceColor', colours{1}, 'MarkerEdgeColor', colours{1}, 'MarkerSize', 4)
errorbar (mdiff, err, 'Color', 'Red', 'LineStyle', 'None')
hold off
end
% create the video writer with 1 fps
writerObj = VideoWriter ( 'AtomCloudShift.avi' );
writerObj.FrameRate = 2;
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open (writerObj);
% write the frames to the video
for i = 1: length (F)
% convert the image to a frame
frame = F (i);    
writeVideo (writerObj, frame);
end
% close the writer object
close (writerObj);
