%%
corrector = ImageCorrection.getInstance();

[image, fLog, filter, filteredImage] = corrector.applyFourierMask('2020-06-04T005311_seq_scanRamseyFringeForHDT2TrapDepthMeasurement.mat', 'UseFFT', true);

FirstImage = squeeze(squeeze(corrector.AverageImages(1,1,:,:)));
%%

figure(1)
colormap(parula)
subplot(2,3,1),imagesc(image); title('Original Image')
colorbar
subplot(2,3,2),imagesc(fLog); title('Fourier Image')
colorbar
subplot(2,3,3),imagesc(filter); title('Fourier mask')
colorbar
subplot(2,3,4),imagesc(filteredImage); title('Fourier mask Applied')
colorbar
subplot(2,3,5),imagesc(image .* filteredImage); title('Corrected Image')
colorbar
subplot(2,3,6); title('Figure of Merit')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Correction for Etaloning', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 15, ...
    'FontWeight', 'bold')
            
%%

figure(2)
subplot(2,2,1)
s1 = surf(FirstImage,'FaceAlpha',1,'Linestyle','none');
title('Original Image')
s1.EdgeColor = 'none';
subplot(2,2,2)
s2 = surf(FirstImage.*filteredImage,'FaceAlpha',1,'Linestyle','none');
title('Corrected Image')
s2.EdgeColor = 'none';
subplot(2,2,3)
s3 = surf(FirstImage,'FaceAlpha',1,'Linestyle','none');
s3.EdgeColor = 'none';
view([0 0])
subplot(2,2,4)
s4 = surf(FirstImage.*filteredImage,'FaceAlpha',1,'Linestyle','none');
s4.EdgeColor = 'none';
view([0 0])
for kk =1:4
    subplot(2,2,4)
    axis tight
end