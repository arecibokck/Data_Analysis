classdef ImageCorrection < handle
    
    properties
        
        pathToFile
        filename
        imagesLoaded = false;
        measObj
    
        ImageToCorrect
        AverageImages
        AverageBackgroundImage
        AverageBackGroundImageFull
        FourierFilteringType
    end
    
    methods % Lifecycle functions
        function this = ImageCorrection()
            scriptpath = matlab.desktop.editor.getActiveFilename;
            this.pathToFile = scriptpath(1:strfind(scriptpath, [filesep '@' class(this)])-1);
        end
    end
    
    methods % Action handlers
        function loadImages(this)
            if ~this.imagesLoaded
                imagefile = load([this.pathToFile filesep this.filename]);
                this.imagesLoaded = true;
            end
            fname = fieldnames(imagefile);
            this.measObj = getfield(imagefile, fname{1});
            [this.AverageImages, this.AverageBackgroundImage, this.AverageBackGroundImageFull] = this.measObj.getAverageImages;
        end
        function ret = applyBuiltInFFTFilter(~, image)
            % normalize the image and conver to doubleI
            image = double(mat2gray(image));

            % Resize the image
            % I = imresize(I, [256 256]);

            % get the size of the image
            %[rows,cols] = size(image);

            % apply FFT
            f = fftshift(fft2(image));

            % used to plot the image
            fLog = log(abs(f));

            % filter by a range based on fLog

            filter = (fLog < .5*max(fLog(:)) ) | (fLog > 0.88*max(fLog(:)) );

            ret = abs(ifft2(f.*filter));

            figure(1)
            colormap(gray)
            subplot(2,2,1),imagesc(image); title('Original Image')
            colorbar
            subplot(2,2,2),imagesc(fLog); title('Fourier Image')
            colorbar
            subplot(2,2,3),imagesc(filter); title('Mask')
            colorbar
            subplot(2,2,4),imagesc(ret); title('Corrected Image')
            colorbar
            annotation('textbox', [0 0.9 1 0.1], ...
                'String', 'Correction for Etaloning', ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 15, ...
                'FontWeight', 'bold')

            
        end
        function ret = applyintensityGradientCorrection(~, image)
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
        function correctEtaloning(this, filename, varargin)
            
            input = inputParser;
            addRequired(input,'filename',@ischar);
            addParameter(input,'ImageToCorrect','Bkg',@ischar);
            addParameter(input,'FourierFilteringType','BuiltIn',@ischar);
            parse(input, filename, varargin{:});
            
            this.filename = input.Results.filename;
            this.ImageToCorrect = input.Results.ImageToCorrect;
            this.FourierFilteringType = input.Results.FourierFilteringType;
            
            if ~this.imagesLoaded
                this.loadImages;
                this.imagesLoaded = true;
            end
            
            switch this.ImageToCorrect
                case 'Bkg'
                    image = squeeze(squeeze(this.AverageBackgroundImage));
            end
            
            switch this.FourierFilteringType 
                case 'BuiltIn'
                    filteredImage = this.applyBuiltInFFTFilter(image);
                case 'Manual'
                    filteredImage = this.applyManualFFTFilter(image);
            end
            
            FirstImage = squeeze(squeeze(this.AverageImages(1,1,:,:)));
            figure
            subplot(2,2,1)
            s1 = surf(FirstImage,'FaceAlpha',0.5);
            title('Original Image')
            s1.EdgeColor = 'none';
            subplot(2,2,2)
            s2 = surf(FirstImage.*filteredImage,'FaceAlpha',0.5);
            title('Corrected Image')
            s2.EdgeColor = 'none';
            subplot(2,2,3)
            s3 = surf(FirstImage,'FaceAlpha',0.5);
            s3.EdgeColor = 'none';
            view([0 0])
            subplot(2,2,4)
            s4 = surf(FirstImage.*filteredImage,'FaceAlpha',0.5);
            s4.EdgeColor = 'none';
            view([0 0])
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
    end
    
end