%% OPLL ACTIVE FOR BOTH ARMS

%Scanning up in frequency
load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL active for both/2019-12-05T170313_seq_scanLaserFrequency.mat');
meas = meas20191205T170313;

%% OPLL ACTIVE FOR BOTH ARMS

%Scanning down in frequency
load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL active for both/2019-12-05T174019_seq_scanLaserFrequency.mat');
meas = meas20191205T174019;

%% OPLL ACTIVE FOR ONE ARM
%Scanning up in frequency

%HDT3 latched ~1Ghz scan
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL active for one arm/2019-12-05T193228_seq_scanLaserFrequency.mat');
% meas = meas20191205T193228;

%HDT1 latched ~2Ghz scan
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL active for one arm/HDT1 latched/2019-12-05T215540_seq_scanLaserFrequency.mat');
% meas = meas20191205T215540;

%HDT3 latched ~2Ghz scan
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL active for one arm/HDT3 latched/2019-12-05T225438_seq_scanLaserFrequency.mat');
% meas = meas20191205T225438;

%HDT1 latched ~2Ghz scan
load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-07_LaserFrequencyScanforPLD/OPLL active for one arm/HDT1 latched/2019-12-07T190344_seq_scanLaserFrequency.mat');
meas = meas20191207T190344;

%HDT3 ~2Ghz scan
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-07_LaserFrequencyScanforPLD/OPLL active for one arm/HDT3 latched/2019-12-07T195630_seq_scanLaserFrequency.mat');
% meas = meas20191207T195630;

%% OPLL ACTIVE FOR ONE ARM
%Scanning down in frequency

%HDT3 latched ~1Ghz scan
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL active for one arm/2019-12-05T184201_seq_scanLaserFrequency.mat');
% meas = meas20191205T184201;

%HDT1 >2Ghz scan
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL active for one arm/HDT1 latched/2019-12-05T221245_seq_scanLaserFrequency.mat');
% meas = meas20191205T221245;

%HDT3 latched >2Ghz scan
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL active for one arm/HDT3 latched/2019-12-05T231032_seq_scanLaserFrequency.mat');
% meas = meas20191205T231032;

%HDT1 ~2Ghz scan
load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-07_LaserFrequencyScanforPLD/OPLL active for one arm/HDT1 latched/2019-12-07T192707_seq_scanLaserFrequency.mat');
meas = meas20191207T192707;

%HDT3 ~2Ghz scan
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-07_LaserFrequencyScanforPLD/OPLL active for one arm/HDT3 latched/2019-12-07T201726_seq_scanLaserFrequency.mat');
% meas = meas20191207T201726;

%% OPLL INACTIVE FOR BOTH ARMS
%Scanning up in frequency
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL inactive for both arms/2019-12-05T185829_seq_scanLaserFrequency.mat');
% meas = meas20191205T185829;

%Both latched >2Ghz scan
load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL inactive for both arms/2019-12-05T233727_seq_scanLaserFrequency.mat');
meas = meas20191205T233727;

%% OPLL INACTIVE FOR BOTH ARMS
%Scanning down in frequency
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL inactive for both arms/2019-12-05T191411_seq_scanLaserFrequency.mat');
% meas = meas20191205T191411;

%Both latched >2Ghz scan
load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-12-05_LaserFrequencyScanforPLD/OPLL inactive for both arms/2019-12-05T235754_seq_scanLaserFrequency.mat');
meas = meas20191205T235754;

%%
Images = struct([]);
figure('Position', [70, 500, 1800, 400])
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
n = 50;
for i = 1:n
    Images{i} = squeeze(meas.runData.images(i,:,:,:));
end

F(n) = struct('cdata',[],'colormap',[]);
mdiff = [];
err = [];
frame = 1;

diffinpixelsX = zeros(1,n);
diffinpixelsY = zeros(1,n);

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
    
    subplot(1,3,1), imagesc(img1)
    colorbar
    title('First Image')
    
    subplot(1,3,2), imagesc(img2)
    colorbar
    title('Second Image')
    
    subplot(1,3,3)
    plot(Image1Sum);
    hold on 
    plot(Image2Sum);
    
    F (frame) = getframe (gcf);
    drawnow
    frame = frame+1;
    if i<n cla; end
end

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
