figure(1);clf
zoom on; grid on; hold on

plot(FirstImageSumArray(21, :))

% for j = 21
% plot(FirstImageSumArray(j, :))
% ImageNumbers{end+1} = num2str(j);
% end

x_vals = 1:size(FirstImageSumArray,2);

height_env = 0.9;
center_env = 300;
width_env = 80;
contrast =0.7;
period=15.2;
phase_shift = 0.6;
model = height_env*exp(-(x_vals-center_env).^2/(2*width_env^2)) .* ((1+contrast*sin(x_vals/period-phase_shift*2*pi))/2);

plot(model);

%%

figure(1);clf
zoom on; grid on; hold on

plot(FirstImageSumArray(1, :))

% for j = 21
% plot(FirstImageSumArray(j, :))
% ImageNumbers{end+1} = num2str(j);
% end

x_vals = 1:size(FirstImageSumArray,2);

height_env = 0.9;
center_env = 300;
width_env = 80;
contrast =0.7;
period=15.2;
phase_shift = 0.2;
model = height_env*exp(-(x_vals-center_env).^2/(2*width_env^2)) .* ((1+contrast*sin(x_vals/period-phase_shift*2*pi))/2);

plot(model);

%%

figure(1);clf
zoom on; grid on; hold on

plot(FirstImageSum(50, :))

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
