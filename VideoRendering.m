fprintf('Rendering video...\n');
    
filename = strcat('output.avi');

%create a video
workingDir = pwd;
outputVideo = VideoWriter(fullfile(workingDir, filename), 'Uncompressed AVI');
outputVideo.FrameRate = 3.5;
open(outputVideo);

for m = 1:17
    openfig(['figure' num2str(m) '.fig']);
    hold on;
    writeVideo(outputVideo, getframe(gcf));
end                
close(outputVideo);
fprintf('Done.\n');