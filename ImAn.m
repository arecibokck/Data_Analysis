
close all;
% meas = meas20190916T174417;
meas = meas20190919T155409;
ImageArrays = meas.runData.images;
NofImages = size(ImageArrays,1);
FirstImageSum = []; periods = []; err=[]; ImageNumbers = {};
for i = 1:NofImages
%   Image = squeeze(squeeze(flipud(ImageArrays(i, 1, :, :))));
%   ImageSum = sum(Image, 1) / max(sum(Image, 1)); 
    [Image,d] = spdiags(imrotate(squeeze(squeeze(flipud(ImageArrays(i, 1, :, :)))), -90));
    ImageSum = transpose(sum(Image, 2) / max(sum(Image, 2))); 
    Offset = min(ImageSum);
    FirstImageSum = vertcat(FirstImageSum, ImageSum - Offset);
end

%%
figure(1);clf
zoom on; grid on; hold on

plot(FirstImageSum(30, :))

% for j = 21
% plot(FirstImageSumArray(j, :))
% ImageNumbers{end+1} = num2str(j);
% end

x_vals = 1:size(FirstImageSum,2);

height_env = 0.9;
center_env = 340;
width_env = 80;
contrast =0.7;
period=103; %95.5;
phase_shift = 0.7;
model = height_env*exp(-(x_vals-center_env).^2/(2*width_env^2)) .* ((1+contrast*sin(2*pi*x_vals/period-phase_shift*2*pi))/2);

plot(model);

% figure(1)
% zoom on; grid on; hold on
% for j = 1:2:NofImages
%     plot(FirstImageSum(j, :))
%     ImageNumbers{end+1} = num2str(j);
% end
% xtix=get(gca,'xtick')';
% set(gca,'xticklabel',num2str(round(xtix/1.7, -1)))
% title(['Imaging in the Dark']);
% xlabel('Lattice Sites'); 
% ylabel('Normalized Flourescence');

% lg = legend(ImageNumbers);
% title(lg, 'Image Numbers')


