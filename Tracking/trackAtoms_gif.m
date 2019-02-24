function trackAtoms(measData,varargin)

%default values
n = 0;
maxdisp = 30; %an estimate of the maximum distance (in pixels) that a particle would move in a single time interval
peak_find_thres = 200; %Threshold for peak find in an image
pol_thres = 2; %Threshold for selecting clean trajectories 

if (~isempty(varargin))
    switch nargin
        case 2
            n = varargin{1};
        case 3
            n = varargin{1};
            maxdisp = varargin{2};
        case 4    
            n = varargin{1};
            maxdisp = varargin{2};
            peak_find_thres = varargin{3};
        case 5
            n = varargin{1};
            maxdisp = varargin{2};
            peak_find_thres = varargin{3};
            pol_thres = varargin{4};
        otherwise 
            error('Invalid number of optional input arguments');   
    end   
end

if n == 0
    
    n_run = length(measData.runData.images(:,1,1,1)); %Extract number of runs in data set
    l_diffpeaks = zeros(1,n_run); %Initialize an array for collecting number of peaks found in the difference image - essentially, counting the # of atoms

    for k = 1:n_run %Iterate through number of runs
        A = squeeze(measData.runData.images(k,1,:,:)); %First image in each run
        B = squeeze(measData.runData.images(k,end,:,:)); %Last Image in each run
        diff = imabsdiff(A, B); %Difference between the 2 images, the least difference corresponds to atoms having returned to their original positions
        diffpeaks = FastPeakFind(diff, peak_find_thres); %Peak find in difference image
        l_diffpeaks(k) = length(diffpeaks'); %collecting number of peaks(atoms) in all difference images, there being one difference image for each run
    end

    [A, I] = min(l_diffpeaks(l_diffpeaks > 0), [], 2); %Least difference corresponds to atoms returning to their initial positions after transport
    n = I(1); %So, this is the run we are interested in
    
end


all_pos = {};
for i = 1:size(measData.runData.images, 2) %Iterate through number of images in each run
    img = squeeze(measData.runData.images(n,i,:,:)); %Take one image
    p = FastPeakFind(img, peak_find_thres); %Find peak in image
    x = p(1:2:end);
    y = p(2:2:end);
    current_pos = [x,y]; %x, y coordinates of each peak in an image
    all_pos{i} = current_pos; %store all coordinates in a cell array, each index corresponds to one image
end

prev_pos = zeros(size(all_pos{1},1), 3); 

%To create an array with coordinates across all images and corresponding
%time stamps (in terms of frame (image) numbers) for the purpose of
%tracking. The format is as required by the tracking algorithm used
%subsequently

for j = 1:length(all_pos)%iterate over coordinates stored in a cell array with indices referring to frame number
    if ~isempty(all_pos{j})
        pos = horzcat(all_pos{j}, repmat(j, size(all_pos{j},1), 1)); %store coordinates with a third column indicating frame number
        if (j > 1)
            pos = cat(1, prev_pos, pos); %concatenate all coordinates across images, with a third column indicating frame number; obtain one 3 column matrix with all coordinates.
        end
        prev_pos = pos;
    end
end

%Coordinates of atoms are scrambled, the matrix has to be unscrambled to
%extract individual trajectories of atoms. This is achieved with some open
%source MATLAB code by John C. Crocker (University of Pennyslvania)
%and David Grier (Emory University) originally written in IDL http://www.physics.emory.edu/faculty/weeks//idl/

result = zeros(length(pos), 4);
try
    result = track(pos, maxdisp); %send 3 column coordinate matrix and specify maximum displacement of atoms across images to track function
catch
    fprintf('Skipping due to excessive combinitorics...Retry with different parameters or image/dataset.\n');
end

%Returned is a 4 column matrix which is the original 3 column matrix but now 
%with an additional column containing a tag specifying each atom trajectory 
%and the original data rows sorted into those trajectories and are in contiguous 
%blocks, with the time variable a monotonically increasing function inside each block. 

%The following code is to extract trajectories from the 4 column result matrix. Only complete trajectories 
%(atoms tracked throughout all the images) are extracted.
trajectories = {}; 
t = 1; 
x = []; y = [];
for i = 1:max(result, [], 4) %iterate through all identified trajectories, the tag is a number and the maximum is the total number of identified trajectories 
    indr = find(result(:,4) == i); %extract contiguous blocks by identifying those rows with the same tag and finding the indices of these rows
    if (length(indr) == size(measData.runData.images, 2)) %Filters only complete trajectories 
        for j = 1:length(indr)
            x = [x result(indr(j),1)];
            y = [y result(indr(j),2)]; %Extract x and y coordinates from the row     
            trajectories{t} = horzcat(x',y'); %collect them in 2 column matrices and store them in a cell array with each index corresponding to a different complete trajectory
        end
        %Eliminate irregular trajectories that are not close to a regular polygon 
        [cx,cy] = centroid(polyshape({x},{y}));%Find centroid
        distances = sqrt((x-cx).^2 + (y-cy).^2);%Find distances of points from centroid
        if var(distances) < pol_thres %Set threshold for how much variance you would expect, ideally it should be 0
            t = t + 1;
        end
        x = []; y = [];
    end
end   

[a, MSGID] = lastwarn();
warning('off', MSGID);

datasetname = inputname(1);
filename_gif = strcat(datasetname(5:end),'_n_',int2str(n), '_Atom_Transport_Tracked.gif');

pix_v = 1:size(measData.runData.images, 3);
pix_h = 1:size(measData.runData.images, 4);

%pix_v = 40:190;
%pix_h = 200:350;

%pix_v = 80:300; 
%pix_h = 100:320;

diff_v =  0;
diff_h =  0;


img = squeeze(measData.runData.images(n,1,pix_v,pix_h));
imagesc(img);
hold on;
for t = 1:length(trajectories)
    x = trajectories{t}(1:1,1);
    y = trajectories{t}(1:1,2);
    plot(x-diff_h, y-diff_v, 'r--', 'LineWidth', 1);
    %plot(x-diff_h, y-diff_v,'LineStyle','-','Marker','o','MarkerFaceColor',[1,1,1],'MarkerSize',3,'Color',[1.,1,1,0.6],'LineWidth',2)
end

axis square;

xt = get(gca, 'XTick');
yt = get(gca, 'YTick');

%xt = 0:(10/.36):220;
%yt = 0:(10/.36):220;

%xt = 0:(10/.36):150;
%yt = 0:(10/.36):150;

set(gca, 'XTick', xt, 'XTickLabel', floor(xt*.36));
set(gca, 'YTick', yt, 'YTickLabel', floor(yt*.36));
set (gca,'YDir','normal');
set(gca,'FontSize',15);
xlabel 'Distance (\mum)'
ylabel 'Distance (\mum)'
gif(filename_gif,'DelayTime',0.4,'LoopCount',50,'frame',gcf)


if t > 1
    fprintf('%i regular trajectory(ies) extracted. Rendering GIF...\n', t-1);
    
    %Plot tracks of atoms on images and create a GIF
    for m = 2:size(measData.runData.images, 2) 
        img = squeeze(measData.runData.images(n,m, pix_v, pix_h));
        imagesc(img);
        hold on;
        for t = 1:length(trajectories)
            x = trajectories{t}(1:m,1);
            y = trajectories{t}(1:m,2);
            plot(x-diff_h, y-diff_v, 'r--', 'LineWidth', 1);
            %plot(x-diff_h, y-diff_v,'LineStyle','-','Marker','o','MarkerFaceColor',[1,1,1],'MarkerSize',3,'Color',[1.,1,1,0.6],'LineWidth',2)
        end
        gif
    end                
    fprintf('Done.\n');
else
    fprintf('No regular trajectory(ies) found. GIF rendering abandoned. Retry with different parameters or image/dataset.\n');
end