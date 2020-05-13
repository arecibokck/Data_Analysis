%%
angles = [-90 -80 -70 -60 -50 -40 -30 -20 -10 0 +10 +20 +30 +40 +50 +60 +70 +80 +90];
powerIPol_TypeI = [12.8 12.7 12.4 12.7 12.8 12.7 12.7 12.4 12.3 12.5 12.7 12.8 12.7 12.7 12.9 12.7 13.1 12.7 12.7];
powerIIPol_TypeI = [11.2 10.8 9.95 8.25 6.6 4.64 2.78 1.22 0.33 0.0065 0.35 1.49 2.95 4.5 6.67 8.75 10.7 10.9 11.2];
ref_TypeI = [0.12 0.15 0.45 0.9 1.53 2.17 2.63 2.98 3.23 4.19 3.85 3.55 2.8 2.2 1.5 1.03 0.58 0.28 0.18];

%%
% angles = [-90 -80 -70 -60 -50 -40 -30 -20 -10 0 +10 +20 +30 +40 +50 +60 +70 +80 +90];
powerIPol_TypeII = [14.0 13.9 14.0 14.1 13.8 14.1 14.0 14.0 13.7 14.3 14.3 14.2 13.7 14.0 14.2 14.0 13.8 14.1 14.1];
powerIIPol_TypeII = [11.3 10.8 9.58 8.2 6.23 4.57 2.87 1.41 0.232 0.00423 0.341 1.42 2.77 4.45 6.57 8.25 9.83 10.7 10.8];
ref_TypeII = [0.93 0.945 1.05 1.15 1.42 1.7 1.9 2.19 2.31 2.2 2.18 1.98 1.83 1.6 1.395 1.218 1.04 0.96 0.88];
%%
% plot(angles, ref, '--*')
% hold on
% plot(angles, powerIPol, '--*')
% plot(angles, powerIIPol, '--*')
% Reflectivity_TypeI = ref_TypeI./powerIPol_TypeI;
% Reflectivity_TypeII = ref_TypeII./powerIPol_TypeII;

extRatio_TypeI = powerIIPol_TypeI./powerIPol_TypeI;
extRatio_TypeII = powerIIPol_TypeII./powerIPol_TypeII;

% plot(angles, Reflectivity_TypeI, '--*')
% hold on
% plot(angles, Reflectivity_TypeII, '--*')

plot(angles, extRatio_TypeI, '--*')
hold on
plot(angles, extRatio_TypeII, '--*')

xlabel('Analyzer Angle')
% ylabel('Measured Power(mW)')
ylabel('Extinction Ratio')
title('Extinction Ratio of Codixx Plate Polarizers')