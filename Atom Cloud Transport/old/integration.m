%% OPLL ACTIVE FOR BOTH ARMS

%Scanning up in frequency

%% OPLL ACTIVE FOR BOTH ARMS

%Scanning down in frequency

%% OPLL ACTIVE FOR ONE ARM
%Scanning up in frequency


%% OPLL ACTIVE FOR ONE ARM
%Scanning down in frequency


%% OPLL INACTIVE FOR BOTH ARMS
%Scanning up in frequency

%% OPLL INACTIVE FOR BOTH ARMS
%Scanning down in frequency

%%
Images = struct([]);
figure('Position', [70, 70, 1200, 900])
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
n = 50;
for i = 1:n
    Images{i} = squeeze(meas.runData.images(i,:,:,:));
end

F(n) = struct('cdata',[],'colormap',[]);
mdiff = [];
err = [];
frame = 1;

diffinpixels = zeros(1,n);

for i = 1:n
    img = Images{i};
    img1 = squeeze(img(1,:,:));
    img1 = img1 - min (img1 (:));
    img1 = img1 / max (img1 (:));
    
    img2 = squeeze(img(2,:,:));
    img2 = img2 - min (img2 (:));
    img2 = img2 / max (img2 (:));
    
    Image1Sum = sum(spdiags(img1));
    Image2Sum = sum(spdiags(img2));
    
    subplot(2,2,1), imagesc(img1)
    colorbar
    title('First Image')
    
    subplot(2,2,2), imagesc(img2)
    colorbar
    title('Second Image')
    
    subplot(2,2,3)
%     plot(Image1Sum);
%     hold on 
%     plot(Image2Sum);

    [cc,clags] = xcorr(Image1Sum,Image2Sum);
    ccp = plot(clags,cc)
    hold on
    [ac,alags] = xcorr(Image1Sum);
    acp = plot(alags,ac)
    legend ([acp, ccp], { 'First Image Auto-Correlation' , 'Cross-Correlation' }, 'Location', 'northeast');
    title('Spatial Intensity Correlation of image summed along diagonal')
    hold off

    
    diffinpixels(1,i) = max(clags) - max(alags);
    
    
    subplot(2,2,4);
    xlabel('Shift (# of Pixels)');
    ylabel('Occurences');
    title('Shift of atom cloud')
    
    hold on

    F (frame) = getframe (gcf);
    drawnow
    frame = frame+1;
    if i<n cla; end
end
h = histfit(diffinpixels);
h(1).FaceColor = colours{1}; ...[0 0.4470 0.7410]
h(1).FaceAlpha = 0.5;
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% % create the video writer with 1 fps
% writerObj = VideoWriter ( 'IntensityCorrelation.avi' );
% writerObj.FrameRate = 2;
% writerObj.Quality = 100;
% % set the seconds per image
% % open the video writer
% open (writerObj);
% % write the frames to the video
% for i = 1: length (F)
% % convert the image to a frame
% frame = F (i);
% writeVideo (writerObj, frame);
% end
% % close the writer object
% close (writerObj);
