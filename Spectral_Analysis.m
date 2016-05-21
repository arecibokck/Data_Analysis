function [x,y] = Spectral_Analysis(stimulus_file, response_file)

fps = 600;

if (isunix) %# Linux, Mac
    [status, result] = system( ['wc -l ', stimulus_file] );
    numlines = str2num(result);

elseif (ispc) %# Windows
    numlines = str2num( perl('countlines.pl', response_file) );

else
    error('...');

end

x = csvread(stimulus_file);

y = csvread(response_file);

data = iddata(detrend(y), detrend(x), (1/fps));

f = logspace(0,1,100);

G = spa(data,numlines,f);

h = bodeplot(G);
showConfidence(h,1)
figure
h = spectrumplot(G);
showConfidence(h,1)
end
