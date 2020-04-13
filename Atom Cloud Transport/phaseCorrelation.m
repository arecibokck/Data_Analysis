%% Efficient subpixel image registration by cross-correlation.
% Registers two images (2-D rigid translation) within a  fraction
% of a pixel specified by the user. Instead of computing a zero-padded FFT
% (fast Fourier transform), this code uses selective upsampling by a
% matrix-multiply DFT (discrete FT) to dramatically reduce computation time and memory
% without sacrificing accuracy. With this procedure all the image points are used to
% compute the upsampled cross-correlation in a very small neighborhood around its peak. This
% algorithm is referred to as the single-step DFT algorithm in [1].
%
% [1] Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
% "Efficient subpixel image registration algorithms," Opt. Lett. 33,
% 156-158 (2008).
%
% -----------------------------------------------------------------------
%
% Copyright (c) 2016, Manuel Guizar Sicairos, James R. Fienup, University of Rochester
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Rochester nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% --------------------------------------------------------------------------

%% Syntax
% The code receives the FFT of the reference and the shifted images, and an
% (integer) upsampling factor. The code expects FFTs with DC in (1,1) so do not use
% fftshift.
%
%    output = dftregistration(fft2(f),fft2(g),usfac);
%
% The images are registered to within 1/usfac of a pixel.
%
% output(1) is the normalized root-mean-squared error (NRMSE) [1] between f and
% g.
%
% output(2) is the global phase difference between the two images (should be
% zero if images are real-valued and non-negative).
%
% output(3) and output(4) are the row and column shifts between f and g respectively.
%
%    [output Greg] = dftregistration(fft2(f),fft2(g),usfac);
%
% Greg is an optional output, it returns the Fourier transform of the registered version of g,
% where the global phase difference [output(2)] is also compensated.

% This registration technique is well suited to compare images that are captured
% in Fourier domain (i.e. to evaluate an image reconstruction by holography
% or phase retrieval) which are strictly band-limited and exhibit the
% wrap-around effect.
%
% Even though the registration code assumes band-limited images that wrap around, we
% have obtained very good results when applying it to
% band-limited microscope images, and aliased imagery. That is when shifting
% the image brings in new content instead of wrapping it around or when the
% images are not band-limited.

%% Image Registration
% dftregistration.m receives the FT of f and g and the upsampling factor.
% The code expects DC of the FTs at (1,1) so don't use fftshift.
%
% We now use the image registration code to register f and g within 0.01
% pixels by specifying an upsampling parameter of 100

% load('\\rua.iap.uni-bonn.de\users\Karthik\nobackup_LabData\AtomCloudCenterOfMass\2019-11-19T192230_seq_transportAtoms.mat')
% load('/home/karthik/Documents/MATLAB/adwin/trunk/Scripts/+Karthik/2019-11-19T192230_seq_transportAtoms.mat')
meas = meas20191119T192230;
mdiff = [];
err = [];
figure('Position', [70, 500, 1800, 400])
% sgtitle('Efficient subpixel image registration by cross-correlation in the Fourier domain');
colours = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
n = 5;
SortedImages = cell (n, 10);
for i = 1:n
    c = 1;
    for j = 0:n:size(meas.runData.images, 1)-1
        SortedImages{i,c} = squeeze(meas.runData.images(i+j,:,:,:));
        c = c + 1;
    end
end

expshift = [1, 2, 4, 8, 16];
frame = 1;
for i = 1:n
    diffinpixelsX =[];
    diffinpixelsY = [];
    for j = 1 : 10
        
        img = SortedImages{i,j};
        f = squeeze(img(1,:,:));
        f = f - min (f (:));
        f = f / max (f (:));
        
        g = squeeze(img(2,:,:));
        g = g - min (g (:));
        g = g / max (g (:));
        
        usfac = 100;
        [output, Greg] = dftregistration(fft2(f),fft2(g),usfac);
        diffinpixelsX(end+1) = abs(output(4));
        diffinpixelsY(end+1) = abs(output(3));
        
        
        subplot(1,3,1), imagesc(f)
        % hold on
        % plot(0.5*(acxpeak+1),0.5*(acypeak+1), '+', 'Color', [1 1 1], 'MarkerSize', 10, 'LineWidth', 2);
        % hold off
        colorbar
        title('First Image')
        
        subplot(1,3,2), imagesc(g)
        % hold on
        % plot(0.5*(ccxpeak+1),0.5*(ccypeak+1), '+', 'Color', [1 1 1], 'MarkerSize', 10, 'LineWidth', 2);
        % hold off
        colorbar
        title('Second Image')
        
        subplot(1,3,3)
        hold on
        axis image
        xlim([0 18]);
        xticks([0 1 2 4 8 16])
        ylim([0 18]);
%         xtix=get(gca,'xtick')';
%         set(gca,'xticklabel', num2str(xtix/1e+03))
%         
        xlabel('Expected Shift');
        ylabel('Measured shift (# of Lattice Sites)');
        title('Shift of atom cloud')
        F (frame) = getframe (gcf);
        drawnow
        frame = frame+1;
    end
    mdiff(end+1) = mean(diffinpixelsY)/1.6;
    err(end+1) = std(diffinpixelsY)/1.6;
    subplot(1,3,3), plot(expshift(i), mdiff(i), 'o', 'MarkerFaceColor', colours{1}, 'MarkerEdgeColor', colours{1}, 'MarkerSize', 4)
    errorbar (expshift(i), mdiff(i), err(i), 'Color', 'Red', 'LineStyle', 'None')
    hold off
end
% create the video writer with 1 fps
writerObj = VideoWriter ( 'AtomCloudShift.avi' );
writerObj.FrameRate = 2;
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open (writerObj);
% write the frames to the video
for i = 1: length (F)
    % convert the image to a frame
    frame = F (i);
    writeVideo (writerObj, frame);
end
% close the writer object
close (writerObj);