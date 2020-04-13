function [timeVector,phaseVector] = frequencyNoise2(samplingTime,numberOfTimes,amplitudes)
% This function provides a phase noise vector for the phase of a laser
% field given some amplitudes of common noise types.
% 
% call:
% [timeVector,phaseVector] = frequencyNoise2(samplingTime,numberOfTimes,amplitudes)
% 
% ------Output------
% timeVector        vector giving the time signal (equidistant, starting
%                   from samplingTime) 
% phaseVector       vector with the noise phase for the respective point in
%                   time of the timeVector
%
% ------Input-------
% samplingTime      sampling time constant, gives the time step size
% numberOfTimes     integer length of the full noise (and time) vector
% amplitudes        amplitudes of the different noises, currently 3
%                   different noises are implemented, so 3 entries need to
%                   be provided:
%                   amplitudes(1)   parametrizes the gaussian noise of the
%                                   frequency spectral density (fsd)
%                   amplitudes(2)   parametrizes flicker (1/f) noise of the
%                                   fsd
%                   amplitudes(3)   parametrizes random walk (1/f^2) noise
%                                   of the fsd
%                   note that the fsd is proportional to f^2 times the
%                   phase spectral density, which is the magnitude squared
%                   of the fourier transform of the phase

timeVector = (1:numberOfTimes)*samplingTime;

phaseVector = amplitudes(1).*filterGeneratedNoise2(numberOfTimes,1,samplingTime) ...
    + amplitudes(2).* filterGeneratedNoise2(numberOfTimes,1.5,samplingTime) ...
    + amplitudes(3).*filterGeneratedNoise2(numberOfTimes,2,samplingTime);

end


