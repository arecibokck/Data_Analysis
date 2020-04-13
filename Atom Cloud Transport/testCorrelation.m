function testCorrelation()
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
figure('Position', [70, 70, 1200, 900])

[coords, img1] = genAtoms2DImage(100, 5);
img1 = img1 - min (img1 (:));
img1 = img1 / max (img1 (:));

img2 = shiftAtoms2DImage(coords, 0, 5, 5);
img2 = img2 - min (img2 (:));
img2 = img2 / max (img2 (:));

autocorr = normxcorr2_general(img1, img1);
autocorr(1:10,:) = NaN;
autocorr(:,1:10) = NaN;
autocorr(:,end-10:end) = NaN;
autocorr(end-10:end,:) = NaN;
[acypeak,  acxpeak] = find (autocorr == max (autocorr (:)));

crosscorr = normxcorr2_general(img2, img1);
crosscorr(1:10,:) = NaN;
crosscorr(:,1:10) = NaN;
crosscorr(:,end-10:end) = NaN;
crosscorr(end-10:end,:) = NaN;
[ccypeak,  ccxpeak] = find (crosscorr == max (crosscorr (:)));

subplot(2,2,1), imagesc(img1)
axis on;
title('First Image');
colorbar

subplot(2,2,2), imagesc(img2)
axis on;
title('Second Image');
colorbar

subplot(2,2,3)
ac = surf(autocorr, 'FaceColor' , colours{1} , 'FaceAlpha' , 0.9, 'EdgeColor' , 'none' );
hold on
cc = surf(crosscorr, 'FaceColor' , colours{3} , 'FaceAlpha' , 0.9, 'EdgeColor' , 'none' );
legend ([ac, cc], { 'First Image Auto-Correlation' , 'Cross-Correlation' }, 'Location', 'best');
title('Feature matching through Spatial Intensity Correlation')
% shading flat
hold off

diffinpixelsX = abs(ccxpeak - acxpeak);
diffinpixelsY = abs(ccypeak - acypeak);
subplot(2,2,4)
hold on
xlabel('Index');
ylabel('Along direction of shift');
title('Shift of image feature')

subplot(2,2,4), plot(diffinpixelsY, 'o', 'MarkerFaceColor', colours{1}, 'MarkerEdgeColor', colours{1}, 'MarkerSize', 4)

end

function [coordinates, circlePixels] = genAtoms2DImage(nAtoms, radiusAtoms)
[columnsInImage, rowsInImage] = meshgrid(1:489, 1:489);
xvals = randi([100, 400], 1, nAtoms);
yvals = randi([100, 400], 1, nAtoms);
circlePixels = zeros(489,489);
coordinates = zeros(length(xvals), 2);
for i = 1:length(xvals)
    j =randi([1,length(xvals)], 1,1);
    k = randi([1,length(yvals)], 1,1);
    temp = (rowsInImage - yvals(k)).^2 + (columnsInImage - xvals(j)).^2 <= radiusAtoms.^2;
    coordinates(i,1) = xvals(j);
    coordinates(i,2) = yvals(k);
    circlePixels = circlePixels + temp;
end
end

function shiftedPixels = shiftAtoms2DImage(coords, xshift, yshift, radiusAtoms)
[columnsInImage, rowsInImage] = meshgrid(1:489, 1:489);
xvals = coords(:,1) + xshift;
yvals = coords(:,2) + yshift;
shiftedPixels = zeros(489,489);
for i = 1:length(xvals)
    temp = (rowsInImage - yvals(i)).^2 + (columnsInImage - xvals(i)).^2 <= radiusAtoms.^2;
    shiftedPixels = shiftedPixels + temp;
end
end
