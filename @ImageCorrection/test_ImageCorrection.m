%%
clc
corrector = ImageCorrection.getInstance();
                                         
[correctedImages, FourierMask] = corrector.correctEtaloning('CombinedMeasurement.mat', ...
                                                            'CorrectIntensityGradient', true, ...
                                                            'UseFFT', true, ...
                                                            'UseGaussian', true, ...      %UseButterworth', true
                                                            'FilterInnerRadius', 20, ...
                                                            'FilterOuterRadius', 100,...
                                                            'GaussianFilterSigma', 20, ... %'ButterworthFilterOrder', 6, ...
                                                            'SubtractDCOffset', true, ... 
                                                            'DCOffset', 70, ...
                                                            'SaveMask', true, ...
                                                            'CorrectAll', true);
                                                        
% - WITH NO INTENSITY GRADIENT CORRECTION & GAUSSIAN FILTER
clc
corrector.CorrectIntensityGradient = false;
corrector.DCOffset = 110;
corrector.UseGaussianFilter = true; 
corrector.Filter.InnerRadius = 10;
corrector.Filter.OuterRadius = 50;
corrector.GaussianFilterSigma = 10;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage))/9;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);
%%

% - WITH NO INTENSITY GRADIENT CORRECTION, MEASURED DCOFFSET MATRIX & GAUSSIAN FILTER
clc
corrector.CorrectIntensityGradient = false;
Temp_file = matfile([corrector.pathToFile filesep 'DCOffset.mat']);
corrector.DCOffset = Temp_file.imageData; %ReadoutNoise 
corrector.UseGaussianFilter = true; 
corrector.Filter.InnerRadius = 10;
corrector.Filter.OuterRadius = 50;
corrector.GaussianFilterSigma = 10;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage))/9;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);
%%

% - WITH INTENSITY GRADIENT CORRECTION & GAUSSIAN FILTER
clc
corrector.CorrectIntensityGradient = true;
corrector.DCOffset = 110;
corrector.UseGaussianFilter = false; 
corrector.UseButterworthFilter = true;
corrector.Filter.InnerRadius = 10;
corrector.Filter.OuterRadius = 50;
corrector.ButterworthFilterOrder = 6;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage))/9;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);

%%

% - WITH INTENSITY GRADIENT CORRECTION, MEASURED DCOFFSET MATRIX & GAUSSIAN FILTER
clc
corrector.CorrectIntensityGradient = true;
Temp_file = matfile([corrector.pathToFile filesep 'DCOffset.mat']);
corrector.DCOffset = Temp_file.imageData; %ReadoutNoise 
corrector.UseGaussianFilter = true; 
corrector.Filter.InnerRadius = 10;
corrector.Filter.OuterRadius = 50;
corrector.GaussianFilterSigma = 10;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage))/9;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);

%%

% -WITH NO INTENSITY GRADIENT CORRECTION & BUTTERWORTH FILTER
clc
corrector.CorrectIntensityGradient = false;
corrector.DCOffset = 110; %ReadoutNoise 
corrector.UseGaussianFilter = false; 
corrector.UseButterworthFilter = true;
corrector.Filter.InnerRadius = 10;
corrector.Filter.OuterRadius = 50;
corrector.ButterworthFilterOrder = 6;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage))/9;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);

%%
% -WITH NO INTENSITY GRADIENT CORRECTION, MEASURED DCOFFSET MATRIX & BUTTERWORTH FILTER
clc
corrector.CorrectIntensityGradient = false;
Temp_file = matfile([corrector.pathToFile filesep 'DCOffset.mat']);
corrector.DCOffset = Temp_file.imageData; %ReadoutNoise 
corrector.UseGaussianFilter = false; 
corrector.UseButterworthFilter = true;
corrector.Filter.InnerRadius = 10;
corrector.Filter.OuterRadius = 50;
corrector.ButterworthFilterOrder = 6;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage))/9;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);

%%
% -WITH INTENSITY GRADIENT CORRECTION & BUTTERWORTH FILTER
clc
corrector.CorrectIntensityGradient = true;
corrector.DCOffset = 110; %ReadoutNoise 
corrector.UseGaussianFilter = false; 
corrector.UseButterworthFilter = true;
corrector.Filter.InnerRadius = 10;
corrector.Filter.OuterRadius = 50;
corrector.ButterworthFilterOrder = 6;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage))/9;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);

%%
%%
% -WITH INTENSITY GRADIENT CORRECTION, MEASURED DCOFFSET MATRIX & BUTTERWORTH FILTER
clc
corrector.CorrectIntensityGradient = true;
Temp_file = matfile([corrector.pathToFile filesep 'DCOffset.mat']);
corrector.DCOffset = Temp_file.imageData; %ReadoutNoise 
corrector.UseGaussianFilter = false; 
corrector.UseButterworthFilter = true;
corrector.Filter.InnerRadius = 10;
corrector.Filter.OuterRadius = 50;
corrector.ButterworthFilterOrder = 6;
imageData = squeeze(squeeze(corrector.AverageBackgroundImage))/9;
[fLog, filter, EtaloningMask] = corrector.createFourierMask(imageData);

%%                                                        

% FourierMask = corrector.correctEtaloning('2020-06-04T005311_seq_scanRamseyFringeForHDT2TrapDepthMeasurement.mat', ...
%                                          'UseFFT', true, ...
%                                          'UseGaussian', true, ...                
%                                          'FilterInnerRadius', 0, ...
%                                          'FilterOuterRadius', 100,...
%                                          'GaussianFilterSigma', 25, ...
%                                          'SubtractDCOffset', true, ... 
%                                          'DCOffset', 10, ...
%                                          'SaveMask', true, ...
%                                          'CorrectAll', false);

% FourierMask = corrector.correctEtaloning('2020-06-08T175236_seq_scanTest3DCooling.mat', ...
%                                          'UseFFT', true, ...
%                                          'UseGaussian', true, ...      %UseButterworth', true
%                                          'FilterInnerRadius', 20, ...
%                                          'FilterOuterRadius', 100,...
%                                          'GaussianFilterSigma', 20, ... %'ButterworthFilterOrder', 6, ...
%                                          'SubtractDCOffset', true, ... 
%                                          'DCOffset', 70, ...
%                                          'SaveMask', false, ...
%                                          'CorrectAll', false);

%%
corrector.saveFourierMask(FourierMask, 'EtaloningMask')
%%
EtaloningMask = corrector.loadFourierMask('EtaloningMask_2020-06-03.mat');
% EtaloningMask = corrector.loadFourierMask('EtaloningMask_2020-05-19.mat');
%%
% - plotting
figure(1)
clf
colormap(parula)
[xvals,yvals,X,Y,RSquared] = FluoImageAnalysis.ImageAnalysis.getPosition(fLog,'Units','um');
fxvals = [-size(imageData,1)-1:size(imageData,1)-1] .* (1 / (2*size(imageData,1)));
fyvals = [-size(imageData,2)-1:size(imageData,2)-1] .* (1 / (2*size(imageData,2)));

subplot(2,3,1),imagesc(xvals,yvals,imageData); title('Original Image')
colorbar
subplot(2,3,2),imagesc(fxvals,fyvals,fLog); title('Fourier Image')
colorbar
subplot(2,3,3),imagesc(fxvals,fyvals,filter); title('Filter')
colorbar
subplot(2,3,4),imagesc(xvals,yvals,EtaloningMask); title('Etaloning mask')
%caxis([0.4 1.6])
colorbar
CorrectedBackground = imageData .* EtaloningMask;
subplot(2,3,5),imagesc(xvals,yvals,CorrectedBackground); title('Corrected Image')
colorbar
subplot(2,3,6)
histogram(CorrectedBackground(:),'LineStyle','none'); title('Image Histogram')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on, grid minor;
sgtitle('Correction for Etaloning','FontSize', 15)
            

% -calculation
% AverageFirstImage = squeeze(mean(squeeze(corrector.AverageImages(:,1,:,:))))-ones(size(imageData))*corrector.DCOffset;
% AverageFirstImage = squeeze(mean(squeeze(corrector.AverageImages(:,:,:))))/9-ones(size(imageData))*corrector.DCOffset;
AverageFirstImage = squeeze(mean(squeeze(corrector.AverageImages(:,:,:))))/9-corrector.DCOffset;
%RescaleFactor = 1; % Factor of 5 gives a very smooth 
%EtaloningMask = (EtaloningMask-1)*RescaleFactor+1;
CorrectedImage = AverageFirstImage.*EtaloningMask;

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
    if kk>1
    caxis([0 1200])
    end
end