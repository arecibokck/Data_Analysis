%%
clc
corrector = ImageCorrection.getInstance();
corrector.correctEtaloning('2020-06-04T005311_seq_scanRamseyFringeForHDT2TrapDepthMeasurement.mat', 'UseFFT', true);

%%

% -calculations
ReadoutNoise = 70;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage));
imageData = imageData-ones(size(imageData))*ReadoutNoise;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);

% - plotting
figure(1)
clf
colormap(parula)
[xvals,yvals,X,Y,RSquared] = FluoImageAnalysis.ImageAnalysis.getPosition(fLog,'Units','um');
subplot(2,3,1),imagesc(xvals,yvals,imageData); title('Original Image')
colorbar
subplot(2,3,2),imagesc(xvals,yvals,fLog); title('Fourier Image')
colorbar
subplot(2,3,3),imagesc(xvals,yvals,filter); title('Fourier mask')
colorbar
subplot(2,3,4),imagesc(xvals,yvals,EtaloningMask); title('Etaloning mask')
%caxis([0.4 1.6])
colorbar
CorrectedBackground = imageData ./ EtaloningMask;
subplot(2,3,5),imagesc(xvals,yvals,CorrectedBackground); title('Corrected Image')
colorbar
subplot(2,3,6)
histogram(CorrectedBackground(:),'LineStyle','none'); title('Image Histogram')
set(gca, 'YScale', 'log')
sgtitle('Correction for Etaloning','FontSize', 15)
            
%%
% -calculations
AverageFirstImage = squeeze(mean(squeeze(corrector.AverageImages(:,1,:,:))))-ones(size(imageData))*ReadoutNoise;
RescaleFactor = 1; % Factor of 5 gives a very smooth 
EtaloningMask = (EtaloningMask-1)*RescaleFactor+1;
CorrectedImage = AverageFirstImage./EtaloningMask;

% - plotting
figure(2)
clf 
subplot(2,3,1)
imagesc(EtaloningMask);
title('Etaloning Mask')

subplot(2,3,2)
imagesc(AverageFirstImage);
title('Original Image')
%caxis([80 220])

subplot(2,3,3)
imagesc(CorrectedImage);
title('Corrected Image')
%caxis([80 220])

subplot(2,2,3)
histogram(AverageFirstImage(:),'LineStyle','none')
set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
%s1 = surf(AverageFirstImage,'FaceAlpha',1,'Linestyle','none');
%view([75 0])
title('Original Image Histogram')

subplot(2,2,4)
histogram(CorrectedImage(:),'LineStyle','none')
set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
%s2 = surf(CorrectedImage(:),'FaceAlpha',1,'Linestyle','none');
%view([75 0])
title('Corrected Image Histogram')

for kk =1:3
    subplot(2,3,kk)
    axis tight
    colormap jet
    colorbar
end