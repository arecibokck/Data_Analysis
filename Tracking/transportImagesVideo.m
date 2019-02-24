global xoffset
global yoffset
global pixelsPerSite


pixelsPerSite = 1.7;

xoffset = 150;
yoffset = 100;


transportImages = 
transportImages_squeezed = squeeze(transportImages);

for i = (1:size(transportImages_squeezed,1)-0)
    
    figure(1); clf;
    
    imagesc(squeeze(transportImages_squeezed(i,:,:)))
    
    axis image
    
    a_h = gca;
    
    sites = (-410:20:400);
    a_h.XTick = sites*pixelsPerSite;
    a_h.XTickLabel = arrayfun(@(x)sprintf('%d',x),sites-xoffset,'UniformOutput',false);
    
    sites = (-400:20:400);
    a_h.YTick = sites*pixelsPerSite;
    a_h.YTickLabel = arrayfun(@(x)sprintf('%d',x),yoffset-sites,'UniformOutput',false);
    
    xlim([fromSitesToPixelsX(-80),fromSitesToPixelsX(80)])
    ylim(sort([fromSitesToPixelsY(-80),fromSitesToPixelsY(80)]))
    
    xrule = a_h.XAxis;
    yrule = a_h.YAxis;
    
    xrule.FontSize=18;
    yrule.FontSize=18;
    
    xlabel('X (lattice sites)','FontSize',26)
    ylabel('Y (lattice sites)','FontSize',26)
    
    drawOctagon(i,166.5-xoffset,-57+yoffset)
    drawOctagon(i,143-xoffset,-34+yoffset)
    drawOctagon(i,120.5-xoffset,-41+yoffset)
    drawOctagon(i,115.5-xoffset,-93+yoffset)
 
    export_fig(sprintf('img%02d.png',i),'-transparent');
    pause(1.3)
    
    
end

function retVal = fromSitesToPixelsX(x)
global pixelsPerSite
global xoffset

retVal = (x+xoffset)*pixelsPerSite;
end

function retVal = fromSitesToPixelsY(y)
global pixelsPerSite
global yoffset
retVal = -(y-yoffset)*pixelsPerSite;
end

function drawOctagon(n,x0,y0)

coordinates = [
    0,0;
    10,0;
    20,0;
    34,-14;
    34,-24;
    34,-34;
    20,-48;
    10,-48;
    0,-48;
    -14,-34;
    -14,-24;
    -14,-14;
    0,0;
    ];

coordinates = coordinates(1:n,:) + [x0,y0];

if size(coordinates,1) == 1
    
    coordinates = [coordinates;coordinates];
    
end

coordinates = [fromSitesToPixelsX(coordinates(:,1)),fromSitesToPixelsY(coordinates(:,2))];

hold('on');

plot(coordinates(:,1),coordinates(:,2),'LineStyle','-','Marker','o','MarkerFaceColor',[1,1,1],'MarkerSize',5,'Color',[1.,1,1,0.6],'LineWidth',2)


end



