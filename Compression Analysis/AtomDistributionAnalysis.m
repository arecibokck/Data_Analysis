% - Create the Analyzer-Object
if true
    %% - create Analyzer-object
    DataFolder = '/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/Data Analysis/Compression Analysis';
    SaveFolder = 'Results';
    options={};
    %options=[options,{'Parameter','Value'}];
    options=[options,{'DataFolder',DataFolder}];   % Default=userInput
    options=[options,{'SaveFolder', SaveFolder}];   % Default=DataFolder
    options=[options,{'SaveString','ToBeFilled'}];   % Default=DefaultName
    options=[options,{'SavePlot',true}];   % Default=DefaultName
    options=[options,{'ShowPlot', true}];   % Default=DataFolder
    Analyzer = Measurements.AnalyzeMeasurement(options);
    disp('done')
    % - get the FileNames
    Analyzer.getFileNames;
    disp('done')
    %% - Choose Measurement
    k = 1;
    meas = Analyzer.loadMeasurement(k);
    colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
end

%% - 2.) get work with fluo Images    
if true
    %% - Get and correct average images
    [AverageImages,AverageBackgroundImage,AverageBackGroundImageFull] = meas.getAverageImages;
    EtaloningMask = helperFunctions.getEtaloningMask();
    ImageArraySize = size(AverageImages);
    CorrectedAverageImages = AverageImages-AverageBackGroundImageFull;
    UncorrectedAverageImages = CorrectedAverageImages;
    CorrectedAverageImages = bsxfun(@times,CorrectedAverageImages,reshape(EtaloningMask,1,1,ImageArraySize(end-1),ImageArraySize(end)));
    %% - Analyze by plotting an iso-ring around a point
    Ratios_Spread = zeros(1,size(CorrectedAverageImages,1));
    Ratios_PeakIntensities = zeros(1,size(CorrectedAverageImages,1));
    %Ratios_FWHM = zeros(1,size(CorrectedAverageImages,1));
    Point = Measurements.Measurement.getDefaultSettings.Origin;
    options = {};
    options = [options,{'Origin',Point}];
    options = [options,{'Units','um'}];
    [~,~,~,~,RSquared] = FluoImageAnalysis.ImageAnalysis.getPosition(zeros(489),options{:});
    m = size(meas.runData.usedScanParameters.WaitTimeHorizontal,1);
    n = size(meas.runData.usedScanParameters.WaitTimeHorizontal,2);
    spread_errors = zeros(1,m);
    peakint_errors = zeros(1,m);
    Images = cell (m, n);
    for i = 1:m
        c = 1;
        for j = 0:m:size(meas.runData.images,1)-1
            Images{i,c} = squeeze(meas.runData.images(i+j,:,:,:));
            c = c + 1;
        end
    end
    for i = 1:m
        for j = 1 : n
            img = Images{i,j};
            backgroundImage = squeeze(img(3,:,:));
            FirstImage = squeeze(img(1,:,:)) - backgroundImage;
            SecondImage = squeeze(img(2,:,:)) - backgroundImage;
            [Points,Values] = AtomCloudDistribution_Rings(RSquared, FirstImage);
            [fit_params_firstimage, ~] = FlatTopGaussFit(Points,Values,[max(Values);  0;  50; 14]);
            initial_peakint = max(Values);
            [Points,Values] = AtomCloudDistribution_Rings(RSquared, SecondImage);
            [fit_params_secondimage, ~] = GaussFit(Points,Values, [max(Values);  0;  20; 0]);
            final_peakint = max(Values);
            spread = [2*fit_params_firstimage(2) + (2*sqrt(2*log(2)) * fit_params_firstimage(3));2*sqrt(2*log(2)) * fit_params_secondimage(3)];
            Ratios_Spread(j) = spread(2)/spread(1);
            Ratios_PeakIntensities(j) = initial_peakint/final_peakint;
        end
        spread_errors(i) = std(Ratios_Spread(1:n))/sqrt(n);
        peakint_errors(i) = std(Ratios_PeakIntensities(1:n))/sqrt(n);
    end
    figure('Position', [39, 445, 1892, 521])
    win(1) = subplot(1, 3, 1);
    win(2) = subplot(1, 3, 2);
    win(3) = subplot(1, 3, 3);
    set(win,'Nextplot','add')
    for Run = 1:size(CorrectedAverageImages,1)
        FirstImage = squeeze(CorrectedAverageImages(Run,1,:,:));
        SecondImage = squeeze(CorrectedAverageImages(Run,2,:,:));
        [Points,Values] = AtomCloudDistribution_Rings(RSquared, FirstImage);
        initial_peakint = max(Values);
        sgtitle('Distribution around the Center of VDT','FontSize',20)
        plot(win(1), Points,Values,'--*')
        grid(win(1), 'on')
        set(get(win(1),'XLabel'), 'String', 'Distance [um]','FontSize',16);
        set(get(win(1),'YLabel'), 'String', 'Average Fluorescence','FontSize',16);
        hold on
        [fit_params_firstimage, fit_model] = FlatTopGaussFit(Points,Values,[max(Values);  0;  50; 14]);
        plot(win(1), Points, fit_model(fit_params_firstimage,Points), 'Color', colours{5});
        [Points,Values] = AtomCloudDistribution_Rings(RSquared, SecondImage);
        plot(win(1), Points,Values,'--*')
        final_peakint = max(Values);
        [fit_params_secondimage, fit_model] = GaussFit(Points,Values, [max(Values);  0;  20; 0]);
        plot(win(1), Points, fit_model(fit_params_secondimage,Points), 'Color', colours{4});
        hold off
        spread = [2*fit_params_firstimage(2) + (2*sqrt(2*log(2)) * fit_params_firstimage(3));2*sqrt(2*log(2)) * fit_params_secondimage(3)];
        TString1 = 'Initial integrated fluorescence in annulus' ;
        TString2 = ['Initial Spread from Fit = FWHM: ' num2str(spread(1))];
        TString3 = 'Final integrated fluorescence in annulus';
        TString4 = ['Final Spread  from Fit = FWHM: ' num2str(spread(2))];
        legend(win(1), {TString1, TString2, TString3, TString4},'FontSize',8)
        
        Ratios_Spread(Run) = spread(2)/spread(1);
        %         plot(win(2), meas.analysisData.MWSpec.xvals(Run), Ratios_Spread(Run), 'o','MarkerSize',3,...
        %             'MarkerEdgeColor',colours{2},...
        %             'MarkerFaceColor',colours{2});
        %
        %         hold on
        %         errorbar(win(2), meas.analysisData.MWSpec.xvals(Run), Ratios_Spread(Run), spread_errors(Run), 'Color', colours{2}, 'LineStyle', 'None')
        grid(win(2), 'on')
        set(win(2), 'xlim', [0 max(meas.analysisData.MWSpec.xvals)]);
        set(win(2), 'ylim', [0.5 0.8]);
        set(get(win(2),'XLabel'), 'String', 'Wait Time \mus','FontSize',16);
        set(get(win(2),'YLabel'), 'String', 'Ratio','FontSize',16);
        legend(win(2), {'Spread - Final:Initial'}, 'Location', 'northeast','FontSize',10);
        
        Ratios_PeakIntensities(Run) = final_peakint/initial_peakint;
        %         plot(win(3), meas.analysisData.MWSpec.xvals(Run), Ratios_PeakIntensities(Run), '*','MarkerSize',3,...
        %             'MarkerEdgeColor',colours{6},...
        %             'MarkerFaceColor',colours{6});
        %         hold on
        %         errorbar(win(3), meas.analysisData.MWSpec.xvals(Run), Ratios_PeakIntensities(Run), peakint_errors(Run), 'Color', colours{6}, 'LineStyle', 'None')
        grid(win(3), 'on')
        set(win(3), 'xlim', [0 max(meas.analysisData.MWSpec.xvals)]);
        set(win(3), 'ylim', [1 1.3]);
        set(get(win(3),'XLabel'), 'String', 'Wait Time \mus','FontSize',16);
        set(get(win(3),'YLabel'), 'String', 'Ratio','FontSize',16);
        legend(win(3), {'Peak Intensity - Initial:Final'}, 'Location', 'northeast','FontSize',10);
        pause(0.2);
    end
    %% - Analyze strips
    Ratios_Spread = zeros(1,size(CorrectedAverageImages,1));
    Ratios_PeakIntensities = zeros(1,size(CorrectedAverageImages,1));
    %Ratios_FWHM = zeros(1,size(CorrectedAverageImages,1));
    Point = Measurements.Measurement.getDefaultSettings.Origin;
    options = {};
    options = [options,{'Origin',Point}];
    options = [options,{'Units','um'}];
    [xvals,yvals,X,Y,RSquared] = FluoImageAnalysis.ImageAnalysis.getPosition(zeros(489),options{:});
    m = size(meas.runData.usedScanParameters.WaitTimeHorizontal,1);
    n = size(meas.runData.usedScanParameters.WaitTimeHorizontal,2);
    spread_errors = zeros(1,m);
    peakint_errors = zeros(1,m);
    Images = cell (m, n);
    for i = 1:m
        c = 1;
        for j = 0:m:size(meas.runData.images,1)-1
            Images{i,c} = squeeze(meas.runData.images(i+j,:,:,:));
            c = c + 1;
        end
    end
    for i = 1:m
        for j = 1 : n
            img = Images{i,j};
            backgroundImage = squeeze(img(3,:,:));
            FirstImage = squeeze(img(1,:,:)) - backgroundImage;
            SecondImage = squeeze(img(2,:,:)) - backgroundImage;
            [Points,Values] = AtomCloudDistribution_Strips(X,Y,FirstImage);
            [fit_params_firstimage, ~] = GaussFit(Points,Values,[max(Values);  0;  50; 14]);
            initial_peakint = max(Values);
            [Points,Values] = AtomCloudDistribution_Strips(X,Y,SecondImage);
            [fit_params_secondimage, ~] = GaussFit(Points,Values, [max(Values);  0;  20; 0]);
            final_peakint = max(Values);
            spread = [2*fit_params_firstimage(2) + (2*sqrt(2*log(2)) * fit_params_firstimage(3));2*sqrt(2*log(2)) * fit_params_secondimage(3)];
            Ratios_Spread(j) = spread(2)/spread(1);
            Ratios_PeakIntensities(j) = initial_peakint/final_peakint;
        end
        spread_errors(i) = std(Ratios_Spread(1:n))/sqrt(n);
        peakint_errors(i) = std(Ratios_PeakIntensities(1:n))/sqrt(n);
    end
    figure('Position', [39, 445, 1892, 521])
    win(1) = subplot(1, 3, 1);
    win(2) = subplot(1, 3, 2);
    win(3) = subplot(1, 3, 3);
    set(win,'Nextplot','add')
    for Run = 1:size(CorrectedAverageImages,1)
        FirstImage = squeeze(CorrectedAverageImages(Run,1,:,:));
        SecondImage = squeeze(CorrectedAverageImages(Run,2,:,:));
        [Points,Values] = AtomCloudDistribution_Strips(X,Y,FirstImage);
        initial_peakint = max(Values);
        sgtitle('Distribution about HDT13','FontSize',20)
        plot(win(1), Points,Values,'--*')
        grid(win(1), 'on')
        set(get(win(1),'XLabel'), 'String', 'Distance [um]','FontSize',16);
        set(get(win(1),'YLabel'), 'String', 'Average Fluorescence','FontSize',16);
        hold on
%       [fit_params_firstimage, fit_model] = VoigtFit(Points,Values, [0.5, 0.01]);
%       plot(win(1), Points, fit_model(Points, fit_params_firstimage(1), fit_params_firstimage(1)) / fit_model(0,fit_params_firstimage(1), fit_params_firstimage(1))*max(Values) , 'Color', colours{4});
        [fit_params_firstimage, fit_model] = GaussFit(Points,Values, [max(Values);  0;  20; 0]);
        plot(win(1), Points, fit_model(fit_params_firstimage,Points), 'Color', colours{4});
        [Points,Values] = AtomCloudDistribution_Strips(X,Y,SecondImage);
        plot(win(1), Points,Values,'--*')
        final_peakint = max(Values);
        [fit_params_secondimage, fit_model] = GaussFit(Points,Values, [max(Values);  0;  20; 0]);
        plot(win(1), Points, fit_model(fit_params_secondimage,Points), 'Color', colours{5});
        hold off
        spread = [(2*sqrt(2*log(2)) * fit_params_firstimage(1));2*sqrt(2*log(2)) * fit_params_secondimage(3)];
        TString1 = 'Initial integrated fluorescence in strip' ;
        TString2 = ['Initial Spread from Fit = FWHM: ' num2str(spread(1))];
        TString3 = 'Final integrated fluorescence in strip';
        TString4 = ['Final Spread  from Fit = FWHM: ' num2str(spread(2))];
        legend(win(1), {TString1, TString2, TString3, TString4},'FontSize',8)
        Ratios_Spread(Run) = spread(2)/spread(1);
        plot(win(2), meas.analysisData.MWSpec.xvals(Run), Ratios_Spread(Run), 'o','MarkerSize',3,...
            'MarkerEdgeColor',colours{2},...
            'MarkerFaceColor',colours{2});
        hold on
        errorbar(win(2), meas.analysisData.MWSpec.xvals(Run), Ratios_Spread(Run), spread_errors(Run), 'Color', colours{2}, 'LineStyle', 'None')
        grid(win(2), 'on')
        set(win(2), 'xlim', [0 max(meas.analysisData.MWSpec.xvals)]);
        %set(win(2), 'ylim', [0 0.8]);
        set(get(win(2),'XLabel'), 'String', 'Wait Time \mus','FontSize',16);
        set(get(win(2),'YLabel'), 'String', 'Ratio','FontSize',16);
        legend(win(2), {'Spread - Final:Initial'}, 'Location', 'northeast','FontSize',10);
        Ratios_PeakIntensities(Run) = final_peakint/initial_peakint;
        plot(win(3), meas.analysisData.MWSpec.xvals(Run), Ratios_PeakIntensities(Run), '*','MarkerSize',3,...
            'MarkerEdgeColor',colours{6},...
            'MarkerFaceColor',colours{6});
        hold on
        errorbar(win(3), meas.analysisData.MWSpec.xvals(Run), Ratios_PeakIntensities(Run), peakint_errors(Run), 'Color', colours{6}, 'LineStyle', 'None')
        grid(win(3), 'on')
        set(win(3), 'xlim', [0 max(meas.analysisData.MWSpec.xvals)]);
        %set(win(3), 'ylim', [1 1.3]);
        set(get(win(3),'XLabel'), 'String', 'Wait Time \mus','FontSize',16);
        set(get(win(3),'YLabel'), 'String', 'Ratio','FontSize',16);
        legend(win(3), {'Peak Intensity - Final:Initial'}, 'Location', 'northeast','FontSize',10);
        pause(0.2);
    end
%%  
%     plot(meas.runData.usedScanParameters.WaitTimeHorizontal(:,1), Ratios_Spread, '--', 'Color', colours{2}, 'HandleVisibility','off');
%     plot(meas.runData.usedScanParameters.WaitTimeHorizontal(:,1), Ratios_PeakIntensities, '--', 'Color', colours{6},'HandleVisibility','off');
%
%     N=300;
%     [fit_params, fit_model]  = DampedCosFit(meas.runData.usedScanParameters.WaitTimeHorizontal(:,1), ...
%                                             Ratios_Spread, ...
%                                             'ParamGuess',  [0.2; 0.0008; 0.2;  11; 0.88], ...
%                                             'LowerBounds', [0; 0.00079; 0.198; 10.9; 0], ...
%                                             'UpperBounds', [max(Ratios_Spread); 0.00081; 0.205; 11.05; 1]);
%     plot(win(2), linspace(0,max(meas.runData.usedScanParameters.WaitTimeHorizontal(:,1)), N), fit_model(fit_params,linspace(0,max(meas.runData.usedScanParameters.WaitTimeHorizontal(:,1)), N)), '--', 'Color', colours{2});
%     old_legend=get(win(2), 'Legend');
%     old_legend.String{2} = 'Fit';
%     legend(win(2),old_legend.String)
%     M=300;
%     [fit_params, fit_model]  = DampedCosFit(meas.runData.usedScanParameters.WaitTimeHorizontal(:,1), ...
%                                             Ratios_PeakIntensities, ...
%                                             'ParamGuess',  [0.05; 0.05; 0.001;  0.004; 0.88], ...
%                                             'LowerBounds', [0; 0.00079; 0.198; 10.9; 0], ...
%                                             'UpperBounds', [max(Ratios_Spread); 0.00081; 0.205; 11.05; 1]);
%     plot(win(3), linspace(0,max(meas.runData.usedScanParameters.WaitTimeHorizontal(:,1)), M), fit_model(fit_params,linspace(0,max(meas.runData.usedScanParameters.WaitTimeHorizontal(:,1)), M)),'--', 'Color', colours{6});
%     old_legend=get(win(3), 'Legend');
%     old_legend.String{2} = 'Fit';
%     legend(win(3),old_legend.String)
    compressionRatios = analysis_data(...
        'alpha',1-0.68,...
        'PlotVerticalBars','no');
    compressionRatios.xvals  = meas.runData.usedScanParameters.WaitTimeHorizontal(:,1); %meas.analysisData.MWSpec.xvals
    compressionRatios.yvals  = Ratios_Spread;
    compressionRatios.stderr = spread_errors;
    compressionRatios.CI    =  compressionRatios.StdDev2CI;
    CompressionRatioFit = FitDataGauss(...
        compressionRatios.xvals(1:end),...
        compressionRatios.yvals(1:end),...
        @(p, x) p(1) .* exp(-p(2)*x) .* cos(2*pi*p(3).*x*1e-3 + p(4))  + p(5),...    % Fit model
        [0.15; 0.001;  1;  0; 0.6],...                 % Param guess
        'alpha',1-0.68,...
        'ParametersDescription',{'Amplitude','Decay Constant', 'Frequency(kHz)', 'Phase', 'Offset'},...
        'FixedParams', [nan,nan,nan,nan,nan], ...
        'FitToolbox','optimization',...
        'LowerBoundParams',[0,0,0,-inf,0],...
        'UpperBoundParams',[inf,inf,inf,inf,1],...
        'StandardErrors',compressionRatios.stderr);
    figure(2);
    subplot(2,1,1);
    clf;
    CompressionRatioFit.doFit;
    compressionRatios.plotData;
    CompressionRatioFit.plotFitModel;
    CompressionRatioFit.printFitReport;
%%
    compressionRatios.xvals = meas.runData.usedScanParameters.WaitTimeHorizontal(:,1);%meas.analysisData.MWSpec.xvals
    compressionRatios.yvals = Ratios_PeakIntensities;
    compressionRatios.stderr = peakint_errors;
    compressionRatios.CI    =  compressionRatios.StdDev2CI;
    CompressionRatioFit = FitDataGauss(...
        compressionRatios.xvals(1:end),...
        compressionRatios.yvals(1:end),...
        @(p, x) p(1) .* exp(-p(2)*x) .* cos(2*pi*p(3).*x*1e-3 + p(4))  + p(5),...    % Fit model
        [0.1; 0.0001; 1;  0; 1.14],...                 % Param guess
        'alpha',1-0.68,...
        'ParametersDescription',{'Amplitude','Decay Constant', 'Frequency', 'Phase', 'Offset'},...
        'FixedParams', [nan,nan,nan,pi,nan], ...
        'FitToolbox','optimization',...
        'LowerBoundParams',[0,0,0,-inf,0],...
        'UpperBoundParams',[inf,inf,inf,inf,inf],...
        'StandardErrors',compressionRatios.stderr);
    
    figure(2);
    subplot(2,1,1);
    clf;
    CompressionRatioFit.doFit;
    compressionRatios.plotData;
    CompressionRatioFit.plotFitModel;
    CompressionRatioFit.printFitReport;
%%  Estimating atom count from fluorescence
    cf = 3.8785e-05;
    m = size(meas.runData.usedScanParameters.WaitTimeHorizontal,1);
    n = size(meas.runData.usedScanParameters.WaitTimeHorizontal,2);
    AtomCounts_FirstImage = zeros(n,1);
    AvgAtomCounts_FirstImage = zeros(m,1);
    FirstImage_count_errors = zeros(m,1);
    AtomCounts_SecondImage = zeros(n,1);
    AvgAtomCounts_SecondImage = zeros(m,1);
    SecondImage_count_errors = zeros(m,1);
    for i = 1:m
        for j = 1 : n
            img = Images{i,j};
            backgroundImage = squeeze(img(3,:,:));
            FirstImage = squeeze(img(1,:,:)) - backgroundImage;
            SecondImage = squeeze(img(2,:,:)) - backgroundImage;
            %[Points,Values] = AtomCloudDistribution(RSquared, FirstImage);
            AtomCounts_FirstImage(j) = sum(sum(FirstImage)) * cf;
            %[Points,Values] = AtomCloudDistribution(RSquared, SecondImage);
            AtomCounts_SecondImage(j) = sum(sum(SecondImage)) * cf;
        end
        AvgAtomCounts_FirstImage(i) = sum(AtomCounts_FirstImage) / n;
        FirstImage_count_errors(i) = std(AtomCounts_FirstImage)/sqrt(n);
        AvgAtomCounts_SecondImage(i) = sum(AtomCounts_SecondImage) / n;
        SecondImage_count_errors(i) = std(AtomCounts_SecondImage)/sqrt(n);
    end
    set(gcf, 'defaultAxesColorOrder', [[0 0 0];[0 0 0]]);
    HoldTimes = meas.runData.usedScanParameters.WaitTimeHorizontal(:,1);
    plot(HoldTimes,AvgAtomCounts_FirstImage, 'o', ...
         'MarkerEdgeColor', colours{1}, ...
         'MarkerFaceColor', colours{1});
    hold on
    errorbar(HoldTimes,AvgAtomCounts_FirstImage, FirstImage_count_errors, 'Color', colours{1}, 'LineStyle', 'none', 'HandleVisibility', 'off');
    plot(HoldTimes,AvgAtomCounts_SecondImage, 'o', ...
        'MarkerEdgeColor', colours{2}, ...
        'MarkerFaceColor', colours{2});
    errorbar(HoldTimes,AvgAtomCounts_SecondImage, SecondImage_count_errors, 'Color', colours{2}, 'LineStyle', 'none', 'HandleVisibility', 'off');
    title('Estimation from integrated fluorescence','FontSize', 14)
    xlabel('Hold Time','FontSize', 14)
    ylabel('Atom Count','FontSize', 14)
    yyaxis right
    survival = AvgAtomCounts_SecondImage./AvgAtomCounts_FirstImage;
    survivalerrors = survival .* sqrt((FirstImage_count_errors ./ AvgAtomCounts_FirstImage).^2 + (SecondImage_count_errors ./ AvgAtomCounts_SecondImage).^2);
    plot(HoldTimes, survival, 'o', ...
        'MarkerEdgeColor', colours{5}, ...
        'MarkerFaceColor', colours{5});
    errorbar(HoldTimes,survival, survivalerrors, 'Color', colours{5}, 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylabel('Survival','FontSize', 14)
    ylim([0.5 0.7])
    legend({'Pre-Compression', 'Post-Compression', 'Ratio'},'FontSize', 14)    
end
%% - helperFunctions
function [fit_params, fit_model]  = DampedCosFit(x,y, varargin)
    
    p = inputParser;
    
    addRequired(p,  'Xvals', @isnumeric);
    addRequired(p,  'Yvals', @isnumeric);
    addParameter(p, 'ParamGuess' , [max(y); 0.001; 2;  0; 0.8],  @isnumeric);
    addParameter(p, 'LowerBounds', [0; 0; 0; 0; 0],@isnumeric);
    addParameter(p, 'UpperBounds', [max(y); inf; inf; inf; 1],@isnumeric);
    parse(p, x, y, varargin{:})
    x = p.Results.Xvals;
    y = p.Results.Yvals;
    fit_guess = p.Results.ParamGuess; 
    lb = p.Results.LowerBounds;
    ub = p.Results.UpperBounds;
    fit_model = @(p, x) p(1) .* exp(-p(2)*x) .* cos(2*pi*p(3).*x + p(4))  + p(5) ; % Function to fit
    fcn = @(p) sum((fit_model(p,x) - y').^2);                                      % Least-Squares cost function
    %fit_params = fminsearch(fcn, fit_guess);                                      % Minimise Least-Squares    
    fit_params = fmincon(fcn, fit_guess,[],[],[],[],lb,ub);
    RsqStatistic = 1 - (sum((y' - fit_model(fit_params,x)).^2) / sum((y' - mean(y)).^2));
    fit_params = [fit_params;RsqStatistic];
end
function [fit_params, fit_model]  = GaussFit(x,y, varargin)
    narginchk(2,3);
    fit_guess = [max(y);  0;  30;  0];
    if nargin==3
       fit_guess = varargin{1};   
    end
    fit_model = @(p, x) p(1) .* exp(-(x-p(2)).^2./(2*p(3)^2)) + p(4); % Function to fit
    fcn = @(p) sum((fit_model(p,x)' - y').^2);                        % Least-Squares cost function
    fit_params = fminsearch(fcn, fit_guess);                          % Minimise Least-Squares
    RsqStatistic = 1 - (sum((y' - fit_model(fit_params,x)').^2) / sum((y' - mean(y)).^2));
    fit_params = [fit_params;RsqStatistic];
end
function [fit_params, fit_model]  = FlatTopGaussFit(x,y, varargin)
    narginchk(2,3);
    fit_guess = [max(y);  0;  30;  0];
    if nargin==3
       fit_guess = varargin{1};   
    end
    fit_model = @(p, x) p(1) .* exp(-(max(x-p(2),0)).^2./(2*p(3)^2)) + p(4); % Function to fit
    fcn = @(p) sum((fit_model(p,x)' - y').^2);                        % Least-Squares cost function
    fit_params = fminsearch(fcn, fit_guess);                          % Minimise Least-Squares
    RsqStatistic = 1 - (sum((y' - fit_model(fit_params,x)').^2) / sum((y' - mean(y)).^2));
    fit_params = [fit_params;RsqStatistic];
end
function [fit_params, fit_model] = VoigtFit(x,y,varargin)
    narginchk(2,3);
    fit_guess = [2, 2];
    if nargin==3
       fit_guess = varargin{1};   
    end
    voigt = @(pos,sigma,gamma) real(fadf((pos+1i*gamma)/(sqrt(2)*sigma)));
    p_fit = fminunc(@(p) sum((y-voigt(x,p(1),p(2))/voigt(0,p(1),p(2))*max(y)).^2),fit_guess,optimoptions(@fminunc,'Display','none'));
    sigma_nu = p_fit(1);
    gamma = p_fit(2);
    fprintf(['The Gaussian contribution has a spectral width:\t\tsigma_nu = %1.3f\n' ...
        'The Lorentzian contribution has a spectral width:\tgamma = %1.3f\n'],p_fit(1),p_fit(2));
    fit_params = [sigma_nu gamma];
    fit_model = voigt;
end
function [Points, Values, varargout] = AtomCloudDistribution_Rings(RSquared, Image)
    Limits = union(linspace(0,20,10),linspace(20,200,50));
    Points = diff(Limits)/2+ Limits(1:end-1);
    Values = zeros(1,numel(Limits)-1);
    for kk = 1:(numel(Limits)-1)
        if (Limits(kk)+Limits(kk+1)) < (size(Image,1)/2)
            Mask = getDistanceMask(sqrt(RSquared),Limits(kk),Limits(kk+1));
            Values(kk) = sum(sum(Image.*Mask))/sum(Mask(:));
        end
    end
    Values = Values(Values~=0);
    Points = Points(1:length(Values));
    yvals = sum(Image.*Mask);
    [max_y_val, idx_of_max] = max(yvals);
    idxR=find(yvals(idx_of_max:end)<=0.5*max_y_val,1,'first')+idx_of_max-1;
    if(isempty(idxR))
        idxR = length(yvals);
    end
    idxL=find(yvals(1:idx_of_max)<=0.5*max_y_val,1,'last');
    if(isempty(idxL))
        idxL = 1;
    end
    FWHM = ((idxR - idxL)/1.69) * (0.866/2);
    if nargout > 2
        varargout{1} = FWHM;
    end
end
function [Points, Values, varargout] = AtomCloudDistribution_Strips(X,Y,Image)
    A = [10,-10];
    B = [20,-20];
    [Distances] = getDistance(X,Y,A,B);
    DistLimits = [10,20];
    Mask = getDistanceMask(Distances,DistLimits(1),DistLimits(2));
    Limits =union(linspace(0,20,10),linspace(20,100,50));
    Limits = union(Limits,-Limits);
    Points = diff(Limits)/2+ Limits(1:end-1);
    Values = zeros(1,numel(Limits)-1);
    for kk = 1:(numel(Limits)-1)
        Mask = getDistanceMask(Distances,Limits(kk),Limits(kk+1));
        Values(kk) = sum(sum(Image.*Mask))/sum(Mask(:));
    end
    yvals = sum(Image.*Mask);
    [max_y_val, idx_of_max] = max(yvals);
    idxR=find(yvals(idx_of_max:end)<=0.5*max_y_val,1,'first')+idx_of_max-1;
    if(isempty(idxR))
        idxR = length(yvals);
    end
    idxL=find(yvals(1:idx_of_max)<=0.5*max_y_val,1,'last');
    if(isempty(idxL))
        idxL = 1;
    end
    FWHM = ((idxR - idxL)/1.69) * (0.866/2);
    if nargout > 2
        varargout{1} = FWHM;
    end
end
function Mask = getDistanceMask(Distances,minDist,maxDist)
Mask = ((Distances >=minDist) & (Distances<maxDist))+0;
end
function [Distances] = getDistance(XGrid,YGrid,Point1,Point2)
% [Distances] = getDistance(XGrid,YGrid,Point1,Point2)
% computes either the Euclidian distance from Point1 (if only one point
% is given) for each point defined by XGrid,YGrid.
% Or the distance from the line through Point1 and Point2
%
% In the latter case the distance has a sign indicating on which side
% of the line a point is lying


% - inputHandling
assert(isequal(size(XGrid),size(YGrid)),'Size of XGrid and YGrid must be equal')
assert(numel(Point1)==2, 'Point1 must be a 2d-vector');

if nargin ==3
    Distances = sqrt((XGrid-Point1(2)).^2+(YGrid-Point1(1)).^2);
    
elseif nargin ==4
    % - input handling
    assert(numel(Point2)==2, 'Point2 must be a 2d-vector');
    A = Point1(:)';
    B = Point2(:)';
    
    NormalVector = ([0,-1;1,0]*(A-B)')'; NormalVector = NormalVector/norm(NormalVector); NormalVektorPoint = NormalVector;
    C = zeros(489,489,2);
    C(:,:,1) = YGrid;
    C(:,:,2) = XGrid;
    
    % - "repmat" of A and B
    ATemp = C*0; ATemp(:,:,1) = ATemp(:,:,1)+A(1); ATemp(:,:,2) = ATemp(:,:,2)+A(2); A = ATemp;clear ATemp;
    BTemp = C*0; BTemp(:,:,1) = BTemp(:,:,1)+B(1); BTemp(:,:,2) = BTemp(:,:,2)+B(2); B = BTemp;clear BTemp;
    NTemp = C*0; NTemp(:,:,1) = NTemp(:,:,1)+NormalVector(1); NTemp(:,:,2) = NTemp(:,:,2)+NormalVector(2); N = NTemp;clear NTemp;
    
    BA = (B-A);
    CA = (C-A);
    
    
    D  = (CA - repmat(sum(CA.*BA,3)./sum(BA.*BA,3),1,1,2).*BA);
    Distances = sum(D.*N ,3);
end
end
function FF = fadf(z)
%     This program file computes the complex error function, also known as
% the Faddeeva function. The algorithmic implementation utilizes the
% approximations based on the Fourier expansion [1, 2] and the Laplace
% continued fraction [3]. The code covers with high-accuracy the entire
% complex plain required for practical applications (see also optimized C++
% source code from the RooFit package in the work [4]).
%     The code remains highly accurate even at vanishing imaginary argument
% y -> 0, where y = Im[z]. The worst detected accuracy is ~10^-13.
%
% REFERENCES
% [1] S. M. Abrarov and B. M. Quine, Efficient algorithmic implementation
%     of the Voigt/complex error function based on exponential series
%     approximation, Appl. Math. Comput., 218 (2011) 1894-1902.
%     http://doi.org/10.1016/j.amc.2011.06.072
%
% [2] S. M. Abrarov and B. M. Quine, On the Fourier expansion method
%     for highly accurate computation of the Voigt/complex error function
%     in a rapid algorithm, arXiv:1205.1768v1 (2012).
%     http://arxiv.org/abs/1205.1768
%
% [3] W. Gautschi, Efficient computation of the complex error function,
%     SIAM J. Numer. Anal., 7 (1970) 187-198.
%     http://www.jstor.org/stable/2949591
%
% [4] T. M. Karbach, G. Raven and M. Schiller, Decay time integrals in
%     neutral meson mixing and their efficient evaluation,
%     arXiv:1407.0748v1 (2014).
%     http://arxiv.org/abs/1407.0748
%
%     The code is written by Sanjar M. Abrarov and Brendan M. Quine, York
% University, Canada, September 2014. Last modifications to the code were
% made on July 2016 (see the file 'readme.txt' for more information).
ind_neg = imag(z)<0; % if some imag(z) values are negative, then ...
z(ind_neg) = conj(z(ind_neg)); % ... bring them to the upper-half plane
FF = zeros(size(z)); % define array
ind_ext  = abs(z)>8; % external indices
ind_band = ~ind_ext & imag(z)<5*10^-3; % narrow band indices
FF(~ind_ext & ~ind_band) = ...
    fexp(z(~ind_ext & ~ind_band)); % internal area. This area is ...
    % ... the most difficult for accurate computation
FF(ind_ext) = contfr(z(ind_ext)); % external area
FF(ind_band) = smallim(z(ind_band)); % narrow band
    function FE = fexp(z,tauM,maxN) % Fourier expansion approximation, ...
        % ... see [1, 2] for more information

        if nargin == 1 % assign default paramenetrs tauM and maxN
            tauM = 12; % the margin value
            maxN = 23; % the upper limit summation integer
        end

        n = 1:maxN; % initiate an array
        aN = 2*sqrt(pi)/tauM*exp(-n.^2*pi^2/tauM^2); % Fourier coefficients

        z1 = exp(1i*tauM*z); % define first repeating array
        z2 = tauM^2*z.^2; % define second repeating array

        FE = sqrt(pi)/tauM*(1 - z1)./z2; % initiate array FE
        for n = 1:maxN
            FE = FE + (aN(n)*((-1)^n*z1 - 1)./(n^2*pi^2 - z2));
        end
        FE = 1i*tauM^2*z/sqrt(pi).*FE;
    end
    function CF = contfr(z) % the Laplace continued fraction approximation

        bN = 11; % initial integer
        bN = 1:bN;
        bN = bN/2;

        CF = bN(end)./z; % start computing from the last bN
        for k = 1:length(bN) - 1
            CF = bN(end-k)./(z - CF);
        end
        CF = 1i/sqrt(pi)./(z - CF);
    end
    function SIm = smallim(z) % approximation at small imag(z)

        ind_0 = abs(real(z))<5*1e-3;

        % If |Re[z]| < 5*1e-3, then:
        SIm(ind_0) = small_z(z(ind_0));

        x = real(z); % define the repeating array
        ind_poles = false(size(x)); % initiate the array of indices

        k = 1; % the counter
        while k <= 23
            % These indices are to avoid the poles that can strongly ...
            % ... deteriorate the accuracy in computation
            ind_poles = ind_poles | abs(x - k*pi/12)<1e-4;
            k = k + 1; % just to increment the counter
        end

        % Else if |Re[z]| >= 5*1e-3, then:
        SIm(~ind_0 & ~ind_poles) = narr_band(z(~ind_0 & ~ind_poles));
        % -----------------------------------------------------------------
        % Note that the margin value tauM in the line below is taken as ...
        % 12.1 instead of the default value 12. This excludes all poles ...
        % in computation even if Im[z] -> 0.
        SIm(~ind_0 & ind_poles) = narr_band(z(~ind_0 & ind_poles),12.1,23);
        % -----------------------------------------------------------------

        function SZ = small_z(z)

            % This equation improves accuracy near the origin. It is ...
            % obtained by the Maclaurin series expansion.

            % Define the repeating arrays
            zP2=z.^2;
            zP4=zP2.^2;
            zP6=zP2.*zP4;

            SZ = (((6 - 6*zP2 + 3*zP4 - zP6).*(15*sqrt(pi) + ...
                1i*z.*(30 + 10*zP2 + 3*zP4)))/(90*sqrt(pi)));
        end

        function NB = narr_band(z,tauM,maxN) % the narrow band

            % This is just an alternative representation of the ...
            % equation (14) from [1].

            if nargin == 1 % define default parameters
                    tauM = 12; % the margin value
                    maxN = 23; % the upper limit summation integer
            end

            n = 1:maxN; % initiate an array
            aN = 2*sqrt(pi)/tauM*exp(-n.^2*pi^2/tauM^2); % The Fourier ...
                                         % ... expansion coeffisients

            z1 = cos(tauM*z); % define first repeating array
            z2 = tauM^2*z.^2; % define second repeating array

            NB = 0; % initiate the array NB
            for n = 1:maxN
                NB = NB + (aN(n)*((-1)^n*z1 - 1)./(n^2*pi^2 - z2));
            end
            NB = exp(-z.^2) - 1i*((z1 - 1)./(tauM*z) - ...
                tauM^2*z/sqrt(pi).*NB);
        end
    end
% Convert for negative imag(z) values
FF(ind_neg) = conj(2*exp(-z(ind_neg).^2) - FF(ind_neg));
end
