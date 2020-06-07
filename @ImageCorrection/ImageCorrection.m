classdef ImageCorrection < handle
    
    properties
        
        pathToFile
        filename
        prev_filename
        save_filename
        areImagesLoaded = false;
        measObj
        FigureOfMerit
        AverageImages
        AverageBackgroundImage
        AverageBackGroundImageFull
        UseFFT
        SubtractDCOffset
        SaveMask
        CorrectAll
        DCOffset = 10;
        NotchFilter = struct('InnerRadius', 20, ...
                             'OuterRadius', 100);
        GaussianFilterSigma = 20;
    end
    
    methods % Lifecycle functions
        function this = ImageCorrection()
            scriptpath = matlab.desktop.editor.getActiveFilename;
            this.pathToFile = scriptpath(1:strfind(scriptpath, [filesep '@' class(this)])-1);
        end
    end
    
    methods % Action handlers
        function loadImages(this)
            imagefile = load([this.pathToFile filesep this.filename]);
            fname = fieldnames(imagefile);
            this.measObj = getfield(imagefile, fname{1});
            [this.AverageImages, this.AverageBackgroundImage, this.AverageBackGroundImageFull] = this.measObj.getAverageImages;
        end
        function [fLog, filter, FourierMask] = createFourierMask(this, image)
            if this.UseFFT
                if this.SubtractDCOffset
                    image = image-ones(size(image))*this.DCOffset;
                end
                f = fftshift(fft2(image));
                fLog = log(abs(f));
                [~,~,~,~,RSquared] = FluoImageAnalysis.ImageAnalysis.getPosition(fLog,'Units','pixel');
                filter = ones(size(f));
                filter = filter & (RSquared >= this.NotchFilter.InnerRadius^2); filter = filter & (RSquared <= this.NotchFilter.OuterRadius^2);
                filter =imgaussfilt(double(filter),this.GaussianFilterSigma);
                FourierMask = abs(ifft2(f.*filter));
                FourierMask =FourierMask./mean(FourierMask(:));
                Lower = 0.6;
                Upper = 1.4;
                FourierMask = min(max(FourierMask,Lower),Upper);
            else
                % image = double(mat2gray(image));
                % image = imresize(image, [256 256]);
                M = size(image,1);
                N = size(image,2);
                w_m = exp(-1i*2*pi/M);
                w_n = exp(-1i*2*pi/N);
                [k,m] = meshgrid(1:M,1:M);
                [l,n] = meshgrid(1:N,1:N);
                DFT = ((M*N)^(-1)) * (w_m .^((k-1).*(m-1))) .* (w_n.^((l-1).*(n-1))); % generate Vandermonde DFT matrix
                f = DFT .* image;
                fLog = log(abs(f));
                [~,~,~,~,RSquared] = FluoImageAnalysis.ImageAnalysis.getPosition(fLog,'Units','pixel');
                filter = ones(size(f));
                filter = filter & (RSquared >= this.NotchFilter.InnerRadius^2); filter = filter & (RSquared <= this.NotchFilter.OuterRadius^2);
                filter =imgaussfilt(double(filter),this.GaussianFilterSigma);
                FourierMask = abs(conj(DFT) .* (f.*filter));
                FourierMask =FourierMask./mean(FourierMask(:));
                Lower = 0.6;
                Upper = 1.4;
                FourierMask = min(max(FourierMask,Lower),Upper);
            end
        end
        function saveFourierMask(this, FourierMask, save_filename)
            save_path = [this.pathToFile(1:strfind(this.pathToFile, [filesep 'trunk'])+5) filesep '+calibrationData'];
            save_file = [save_filename '_' strtok(this.measObj.date) '.mat'];
            save([save_path filesep save_file], 'FourierMask');
        end
        function ret = loadFourierMask(this, load_file)
            load_path = [this.pathToFile(1:strfind(this.pathToFile, [filesep 'trunk'])+5) filesep '+calibrationData'];
            Temp_file = matfile([load_path filesep load_file]);
            ret = Temp_file.FourierMask;
        end
        function varargout = correctEtaloning(this, filename, varargin)
            input = inputParser;
            addRequired(input,'filename',@ischar);
            addParameter(input,'UseFFT', true,@islogical);
            addParameter(input,'NotchFilterInnerRadius', this.NotchFilter.InnerRadius,@isnumeric);
            addParameter(input,'NotchFilterOuterRadius', this.NotchFilter.OuterRadius,@isnumeric);
            addParameter(input,'GaussianFilterSigma', this.GaussianFilterSigma,@isnumeric);
            addParameter(input,'SubtractDCOffset', true,@islogical);
            addParameter(input,'DCOffset', this.DCOffset,@isnumeric);
            addParameter(input,'SaveMask', true,@islogical);
            addParameter(input,'CorrectAll', false,@islogical);
            parse(input, filename, varargin{:});
            this.filename = input.Results.filename;
            this.UseFFT = input.Results.UseFFT;
            this.NotchFilter.InnerRadius = input.Results.NotchFilterInnerRadius;
            this.NotchFilter.OuterRadius = input.Results.NotchFilterOuterRadius;
            this.GaussianFilterSigma = input.Results.GaussianFilterSigma;
            this.SubtractDCOffset = input.Results.SubtractDCOffset;
            this.DCOffset = input.Results.DCOffset;
            this.SaveMask = input.Results.SaveMask;
            this.CorrectAll = input.Results.CorrectAll;
            if ~this.areImagesLoaded
                this.loadImages;
                this.areImagesLoaded = true;
            end
            if ~isempty(this.prev_filename) && ~strcmp(this.prev_filename, this.filename)
                this.loadImages;
                this.areImagesLoaded = true;
            end
            
            BkgImage = squeeze(squeeze(this.AverageBackgroundImage));
            
            [~, ~, FourierMask] = this.createFourierMask(BkgImage);
            
            if this.CorrectAll
                if nargout <= 2
                    FirstImages = squeeze(squeeze(this.measObj.runData.images(:,1,:,:)));
                    SecondImages = squeeze(squeeze(this.measObj.runData.images(:,2,:,:)));
                    BackgroundImages = squeeze(squeeze(this.measObj.runData.images(:,3,:,:)));
                    for Index = 1:this.measObj.runData.currentRunIndex
                        FirstImages(Index, :, :) = squeeze(FirstImages(Index, :, :)) ./ FourierMask;
                        SecondImages(Index, :, :) = squeeze(SecondImages(Index, :, :)) ./ FourierMask;
                        BackgroundImages(Index, :, :) = squeeze(BackgroundImages(Index, :, :)) ./ FourierMask;
                    end
                    varargout{1} = struct('CorrectedFirstImages',  FirstImages, ...
                                          'CorrectedSecondImages', SecondImages, ...
                                          'CorrectedBkgImages',  BackgroundImages);
                end
                
                if nargout == 2 
                    varargout{2} = FourierMask;
                end
                
            elseif nargout == 1     
                varargout{1} = FourierMask;
            end
            
            if this.SaveMask
                this.saveFourierMask(FourierMask, 'EtaloningMask');
            end
            
            this.prev_filename = this.filename;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Methods (Static)
    methods (Static)
        % Creates an Instance of Class, ensures singleton behaviour (that there
        % can only be one Instance of this class
        function singleObj = getInstance(varargin)
            % Creates an Instance of Class, ensures singleton behaviour
            persistent localObj;
            if isempty(localObj) || ~isvalid(localObj)
                localObj =  ImageCorrection(varargin{:});
            end
            singleObj = localObj;
        end
        function ret = applyintensityGradientCorrection(image)
            x = [size(image,2)/2 size(image,2)];
            y = [size(image,1) size(image,1)/2];
            z = improfile(image,x,y);

            xp = [0         size(image,2)/2];
            yp = [size(image,1)/2         0];
            zp = improfile(image,xp,yp);

            figure
            clf
            subplot(2,1,1)
            imagesc(image)
            hold on
            plot(x,y,'r')
            plot(xp,yp,'r')
            subplot(2,1,2)
            plot(transpose(floor(size(image,2)/2):size(image,2)),z)
            hold on
            plot(transpose(0:floor(size(image,2)/2)+1),zp)
            y = vertcat(zp, z);
            y = y(~isnan(y));
            x = transpose(linspace(0, length(y), length(y)));
            yu = max(y);
            yl = min(y);
            yr = (yu-yl);                                                    % Range of �y�
            yz = y-yu+(yr/2);
            zx = x(yz .* circshift(yz,[1 0]) <= 0);                          % Find zero-crossings
            per = 2*mean(diff(zx));                                          % Estimate period
            ym = nanmean(y);                                                 % Estimate offset
            model = @(b,x)  b(1).*x + b(2);      % Function to fit
            fcn = @(b) sum((model(b,x) - y).^2);                             % Least-Squares cost function
            s = fminsearch(fcn, [yr;  per;  -1;  ym]);                       % Minimise Least-Squares
            plot(x,model(s,x))

            ret = zeros(size(image,1),size(image,2));
            for row=1:size(ret,1)
                for col=1:size(ret,2)
                    ret(row,col) = ret(row,col) -  s(1)*(row-1);
                end
            end
        end 
    end
end