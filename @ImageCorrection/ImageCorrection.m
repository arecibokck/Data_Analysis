classdef ImageCorrection < handle
    
    properties
        
        pathToFile
        filename
        areImagesLoaded = false;
        measObj
        FigureOfMerit
        AverageImages
        AverageBackgroundImage
        AverageBackGroundImageFull
        UseFFT
    end
    
    methods % Lifecycle functions
        function this = ImageCorrection()
            scriptpath = matlab.desktop.editor.getActiveFilename;
            this.pathToFile = scriptpath(1:strfind(scriptpath, [filesep '@' class(this)])-1);
        end
    end
    
    methods % Action handlers
        function loadImages(this)
            if ~this.areImagesLoaded
                imagefile = load([this.pathToFile filesep this.filename]);
                this.areImagesLoaded = true;
            end
            fname = fieldnames(imagefile);
            this.measObj = getfield(imagefile, fname{1});
            [this.AverageImages, this.AverageBackgroundImage, this.AverageBackGroundImageFull] = this.measObj.getAverageImages;
        end
        function [fLog, filter, FourierMask] = createFourierMask(this, image)
            if this.UseFFT
                % apply FFT
                f = fftshift(fft2(image));
                % used to plot the image
                fLog = log(abs(f));
                %fLog = abs(f)./max(abs(f(:)));
                [~,~,~,~,RSquared] = FluoImageAnalysis.ImageAnalysis.getPosition(fLog,'Units','pixel');
                % filter by a range based on fLog
                filter = ones(size(f));
                %filter = (fLog >= 0.20*max(fLog(:)) ) & (fLog <= 0.30*max(fLog(:)) );
                filter = filter & (RSquared >= 20^2);
                filter = filter & (RSquared <= 100^2);
                filter =imgaussfilt(double(filter),20);
                FourierMask = abs(ifft2(f.*filter));
                %FourierMask =imgaussfilt(FourierMask,5);
                FourierMask =FourierMask./mean(FourierMask(:));
                Lower = 0.6;
                Upper = 1.4;
                FourierMask = min(max(FourierMask,Lower),Upper);
            else
                % normalize the image and conver to doubleI
                image = double(mat2gray(image));
                % Resize the image
                % image = imresize(image, [256 256]);
                % generate Vandermonde DFT matrix
                M = size(image,1);
                N = size(image,2);
                w_m = exp(-1i*2*pi/M);
                w_n = exp(-1i*2*pi/N);
                [k,m] = meshgrid(1:M,1:N);
                [l,n] = meshgrid(1:N,1:N);
                DFT = ((M*N)^(-1)) * (w_m .^((k-1).*(m-1))) .* (w_n.^((l-1).*(n-1)));
                % apply DFT
                f = DFT .* image;
                % used to plot the image
                fLog = log(abs(f));
                % filter by a range based on fLog
                filter = (fLog < .5*max(fLog(:)) ) | (fLog > 0.9*max(fLog(:)) );
                FourierMask = abs(conj(DFT) .* (f.*filter));
            end
        end
        function correctEtaloning(this, filename, varargin)
            
            input = inputParser;
            addRequired(input,'filename',@ischar);
            addParameter(input,'UseFFT', true,@islogical);
            parse(input, filename, varargin{:});
            
            this.filename = input.Results.filename;
            this.UseFFT = input.Results.UseFFT;
            
            if ~this.areImagesLoaded
                this.loadImages;
                this.areImagesLoaded = true;
            end
            
%             BkgImage = squeeze(squeeze(this.AverageBackgroundImage));
%             
%             [~, ~, FourierMask] = this.createFourierMask(BkgImage);
%             
%             FirstImages = squeeze(squeeze(this.measObj.runData.images(:,1,:,:))) ./ FourierMask;
%             
%             SecondImages = squeeze(squeeze(this.measObj.runData.images(:,2,:,:))) ./ FourierMask;
%             
%             BackgroundImages = squeeze(squeeze(this.measObj.runData.images(:,3,:,:))) ./ FourierMask;
%             
%             correctedImages = struct('CorrectedFirstImages',  FirstImages, ...
%                                      'CorrectedSecondImages', SecondImages, ...
%                                      'CorrectedBkgImages',  BackgroundImages);
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