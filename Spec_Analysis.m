function [x,y] = Spec_Analysis(stimulus_file, response_file, folder, name)

fps = 600;

numlines = length(stimulus_file);

x = csvread(stimulus_file);

y = csvread(response_file);

data = iddata(detrend(y), detrend(x), (1/fps));

f = logspace(0,2,100);

G = spa(data,numlines,f);

set(0,'defaultaxesFontName', 'CMU Serif Roman')
set(0,'defaultaxesFontSize', 12)
h = bodeplot(G);
showConfidence(h,1)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 6];
img_name = strcat(strcat(folder,'\'), name);
print ('-dpng', img_name);
pause('on');
pause(1);
close(gcf);

%figure
%h = spectrumplot(G);
%showConfidence(h,1)

end
