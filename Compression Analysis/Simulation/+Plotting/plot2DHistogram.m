function plot2DHistogram(x,y,varargin)
% hist2d(x,y,varargin) this plots the phase-space-density for ... 
    
    p = inputParser;
    addRequired(p,  'Positions',            @isnumeric)
    addRequired(p,  'Velocities',           @isnumeric)
    addParameter(p, 'nbins',           10,  @isscalar)
    addParameter(p, 'PositionLimits',  [-50 50],  @(x) isDataArray(x) && isnumeric(x))
    addParameter(p, 'VelocityLimits',  [-50 50],  @(x) isDataArray(x) && isnumeric(x))
    addParameter(p, 'CountDensity', false,  @islogical)
    parse(p,x,y,varargin{:})
    
    x = p.Results.Positions;
    y = p.Results.Velocities;
    nbins = p.Results.nbins;
    plim  = p.Results.PositionLimits;
    vlim  = p.Results.VelocityLimits;
    cdflag  = p.Results.CountDensity;
    
    xedges = linspace(min(plim),max(plim),nbins+1);
    yedges = linspace(min(vlim),max(vlim),nbins+1);
    
    % computes bins vectors (as middle of each edges couples)
    xbins = mean(cat(1,xedges(1:end-1),xedges(2:end)));
    ybins = mean(cat(1,yedges(1:end-1),yedges(2:end)));

    % computes bins width vectors and area matrix
    xbw = diff(xedges);
    ybw = diff(yedges);
    [xx,yy] = meshgrid(xbw,ybw);
    a = xx.*yy;
    
    % initiate the result matrix
    n = zeros(length(ybins),length(xbins));
    % main loop to fill the matrix with element counts
    for i = 1:size(n,1)
        k = find(y >= yedges(i) & y < yedges(i+1));
        for j = 1:size(n,2)
            n(i,j) = length(find(x(k) >= xedges(j) & x(k) < xedges(j+1)));
        end
    end
    
    % normalize options
    if cdflag
        n = n./a; % count density
    end
    imagesc(xbins,ybins,n)
    %hold on
    %plot(x,y,'.k','MarkerSize',10)
    %hold off	
end