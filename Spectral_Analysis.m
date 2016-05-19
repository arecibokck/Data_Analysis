function [x,y] = Spectral_Analysis(input_file, output_file)

fps = 600;

if (isunix) %# Linux, Mac
    [status, result] = system( ['wc -l ', input_file] );
    numlines = str2num(result);

elseif (ispc) %# Windows
    numlines = str2num( perl('countlines.pl', output_file) );

else
    error('...');

end

x = csvread(input_file);

y = csvread(output_file);

data = iddata(y, x, (1/fps));

disp(data);

f = logspace(0,1,100);

G = spa(data,numlines,f);

h = bodeplot(G);
showConfidence(h,1)
figure
h = spectrumplot(G);
showConfidence(h,1)
end