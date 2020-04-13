function noise = filterGeneratedNoise2(vecSize,exponent,ts)
% This function generates 1/f^exponent noise by spectral filtering of a
% white noise vector. 
%
% call:
% noise = filterGeneratedNoise2(vecSize,exponent)
%
% -------- Output --------
% noise         noise filled vector
% 
% -------- Input ---------
% vecSize       number of elements of the output noise vector
% exponent      exponent specifying the 1/f^order noise

% make Gaussian noise
noise = randn(1,vecSize);

% fourier transform it
noiseFT = fft(noise);
% normalize for proper 1 Hz^2/Hz @ 1 Hz
noiseFT = noiseFT*sqrt(vecSize)/(sqrt(pi)/2);

% apply a filter
VecSize = ceil((numel(noiseFT))/2);

if mod(numel(noiseFT),2)
    freqs = 1/ts * (1:(VecSize-1)) /(VecSize-1);
    noiseFThalf = noiseFT(2:(VecSize))./(freqs.^exponent);

    noiseFTfull = [0 noiseFThalf,fliplr(conj(noiseFThalf(1:end)))];
else
    freqs = 1/ts * (1:VecSize) /VecSize;
    noiseFThalf = noiseFT(2:(VecSize+1))./(freqs.^exponent);
    
    noiseFTfull = [0 noiseFThalf,fliplr(conj(noiseFThalf(1:end-1)))];
end

% back-transform
noiseIFT = ifft(noiseFTfull);

if 1e3*sum(abs(imag(noiseIFT)))>sum(abs(real(noiseIFT)))
    warning('There is something wrong with the symmetry of the fourier transform in the noise generation.');
end

noise = real(noiseIFT);

end
