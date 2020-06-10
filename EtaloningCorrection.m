%%
clc
DoReload = true;
if DoReload
    corrector = ImageCorrection.getInstance();
    
    %Place all data files in same folder as this script
    DarkNoiseBackgroundData = '2020-06-08T175236_seq_scanTest3DCooling';
    EtaloningMaskTemplateData = '2020-06-10T173240_seq_scanRamseyFringeForVDTTrapDepthMeasurement';
    AtomImagesData = 'CombinedMeasurement.mat';
    ImageRescaleFactor = 9;
    
    %Measure and save Dark noise background measurement as DC Offset matrix 
    %meas = load(DarkNoiseBackgroundData);
    %fname = fieldnames(meas);
    %DarkNoiseBackgroundDataObj = getfield(meas, fname{1});
    %AverageImages = DarkNoiseBackgroundDataObj.getAverageImages;
    %AverageFirstImage = squeeze(AverageImages(1,:,:));
    %AverageSecondImage = squeeze(AverageImages(2,:,:));
    %AverageBkgImage = squeeze(AverageImages(3,:,:));
    %AverageAllImages = (AverageFirstImage + AverageSecondImage + AverageBkgImage)/3;
    %save([corrector.pathToFile filesep 'DCOffset.mat'], 'AverageAllImages');
            
    %Dataset of atom images to correct
    meas = load(AtomImagesData);
    fname = fieldnames(meas);
    AtomImageDataObj = getfield(meas, fname{1});
    AverageImages = AtomImageDataObj.getAverageImages/ImageRescaleFactor;
    AverageFirstImage = squeeze(AverageImages(1,:,:));
    AverageSecondImage = squeeze(AverageImages(2,:,:));
    AverageBkgImage = squeeze(AverageImages(3,:,:));
    AverageAllImages = (AverageFirstImage + AverageSecondImage + AverageBkgImage)/3;
    
    corrector.loadImages(EtaloningMaskTemplateData);
end

TypeOfFilter = 'Butterworth';
ImageTemplateToUse = 4; % 1: Avg First Image; 2: Avg Second Image; 3: Avg Bkg Image; 1: Avg All Image
RescaleColorAxis = false;
CustomColorScale = [80,320];
FinalFilterSize = 10;
%DCOffset_file = matfile([corrector.pathToFile filesep 'DCOffset.mat']);
ImageToCorrect = AverageFirstImage;
PlotMaskGenerationFromTemplate = true;
PlotAtomImageCorrection = true;
PlotSurfacePlots = true;

% - Generate mask from image data
switch TypeOfFilter
    case 'Gaussian'
        [fLog, filter, EtaloningMask] = corrector.generateMask(EtaloningMaskTemplateData, ...
                                                               'UseImage'    ,ImageTemplateToUse,  ...
                                                               'CorrectIntensityGradient', false,  ...
                                                               'UseFFT',                    true,  ...
                                                               'UseGaussian',               true,  ...      
                                                               'UseButterworth',           false,  ...
                                                               'FilterInnerRadius',           15,  ...
                                                               'FilterOuterRadius',           30,  ...
                                                               'GaussianFilterSigma',          5,  ... 
                                                               'SubtractDCOffset',          true,  ... 
                                                               'DCOffset',                  3000,  ...
                                                               'SaveMask',                  true);
        
    case 'Butterworth'
        [fLog, filter, EtaloningMask] = corrector.generateMask(EtaloningMaskTemplateData, ...
                                                               'UseImage'    ,ImageTemplateToUse, ...
                                                               'CorrectIntensityGradient', false, ...
                                                               'UseFFT',                    true, ...
                                                               'UseGaussian',              false, ...      
                                                               'UseButterworth',            true, ...
                                                               'FilterInnerRadius',          0.4, ...
                                                               'FilterOuterRadius',           10, ...
                                                               'ButterworthFilterOrder',       3, ...
                                                               'SubtractDCOffset',          true, ... 
                                                               'DCOffset',                    15, ...
                                                               'SaveMask',                  true);        
end

% - Plotting
AvgFirst  = squeeze(corrector.AverageImages(1,:,:)); 
AvgSecond = squeeze(corrector.AverageImages(2,:,:));
AvgBkg    = squeeze(squeeze(corrector.AverageBackgroundImage));
AvgAll    = (AvgFirst + AvgSecond + AvgBkg)/3;

switch ImageTemplateToUse
    case 4
        TemplateImage = AvgAll;
    case 3
        TemplateImage = AvgBkg;
    case 2
        TemplateImage = AvgSecond;
    case 1
        TemplateImage = AvgFirst;
end 

if PlotMaskGenerationFromTemplate
    figure(1)
    clf
    colormap(jet)
    [xvals,yvals,X,Y,RSquared] = FluoImageAnalysis.ImageAnalysis.getPosition(fLog,'Units','um');
    fxvals = [-size(TemplateImage,1)-1:size(TemplateImage,1)-1] .* (1 / (2*size(TemplateImage,1)));
    fyvals = [-size(TemplateImage,2)-1:size(TemplateImage,2)-1] .* (1 / (2*size(TemplateImage,2)));

    subplot(2,3,1),imagesc(xvals,yvals,TemplateImage); title('Original Image')
    colorbar
    subplot(2,3,2),imagesc(fxvals,fyvals,fLog); title('Fourier Image')
    colorbar
    subplot(2,3,3),imagesc(fxvals,fyvals,filter); title('Filter')
    colorbar
    subplot(2,3,4),imagesc(xvals,yvals,EtaloningMask); title('Etaloning mask')
    caxis([0.6 1.4])
    colorbar
    CorrectedBackground = TemplateImage .* EtaloningMask;
    subplot(2,3,5),imagesc(xvals,yvals,CorrectedBackground); title('Corrected Image')
    %caxis([4000 5000])
    colorbar
    subplot(2,3,6)
    histogram(CorrectedBackground(:),'LineStyle','none'); title('Image Histogram')
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    grid on, grid minor;
    sgtitle('Correction for Etaloning','FontSize', 15)
end

if PlotAtomImageCorrection
    % - Atom Image Correction
    CorrectedImage = ImageToCorrect.*EtaloningMask;

    figure(2)
    clf 
    subplot(2,3,1)
    imagesc(EtaloningMask);
    title('Etaloning Mask')

    subplot(2,3,2)
    imagesc(AverageFirstImage);
    title('Original Image')
    if RescaleColorAxis;caxis(CustomColorScale);end

    subplot(2,3,3)
    imagesc(CorrectedImage);
    title('Corrected Image')

    if RescaleColorAxis;caxis(CustomColorScale);end

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
        %caxis([0 1200])
        end
    end
end

if PlotSurfacePlots
    % - Surf plots
    figure(3); 
    clf
    subplot(2,2,1)
    surf(imgaussfilt(TemplateImage,FinalFilterSize),'LineStyle','none');
    % A=5379;B=4300; PeakValley = (A-B)/(A+B)*2;
    % title(sprintf('Uncorrected Background; Flatness = %0.2f', PeakValley))
    title('Uncorrected Background')
    subplot(2,2,2)
    surf(imgaussfilt(CorrectedBackground,FinalFilterSize),'LineStyle','none');
    % A=5004;B=4880; PeakValley = (A-B)/(A+B)*2;
    % title(sprintf('Corrected Background; Flatness = %0.2f', PeakValley))
    title('Corrected Background')
    subplot(2,2,3)
    surf(imgaussfilt(AverageFirstImage,FinalFilterSize),'LineStyle','none');
    % A=1662;B=1326; PeakValley = (A-B)/(A+B)*2;
    % title(sprintf('Uncorrected Atom Image; Flatness = %0.2f', PeakValley))
    title('Uncorrected Atom Image')
    subplot(2,2,4)
    surf(imgaussfilt(CorrectedImage,FinalFilterSize),'LineStyle','none');
    % A=1646;B=1519; PeakValley = (A-B)/(A+B)*2;
    % title(sprintf('Corrected Atom Image; Flatness = %0.2f', PeakValley))
    title('Corrected Atom Image')
    for kk=1:4
        subplot(2,2,kk)
        view([-74 38]);
        colorbar
    end
    sgtitle('Corrected Images','FontSize',20)
end