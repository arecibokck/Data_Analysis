% This function outputs frequency vs. amplitude plot of an input given time series
% The function requires the following inputs
% signal=times series data
% sam=frequency of sampling

function fast_FT(stimulus, response, sam, folder, name)

    %Peak at 0 Hz due to DC mean component
    signal_s = detrend(csvread(stimulus)); %detrend data to remove DC Mean Component (subtract average from each 
    signal_r = detrend(csvread(response)); %data-point) - Now there is no peak at 0 Hz 

    Fs = sam;

    n = length(signal_s);
    L = n; % Window Length of FFT
    nfft = 2^nextpow2(L); % Transform length
    padding = round((nfft - L)/2); %zero padding length on either side of signal

    msignalsdft = Signal_Mod(signal_s, L, nfft, padding); 
    msignalrdft = Signal_Mod(signal_r, L, nfft, padding); 

    f = Fs*(-nfft/2:nfft/2-1)/nfft; %Frequency Range

    dt = 1/Fs; %Frequency resolution


    figure(1)
    set(0,'defaultaxesFontName', 'CMU Serif Roman')
    set(0,'defaultaxesFontSize', 12)
    subplot(2,2,1)
    plot(((0:(n-1)) * dt),signal_s)
    title('Stimulus');
    xlabel('Time (s)'); 
    ylabel('Amplitude (degrees)');
    subplot(2,2,3)  
    plot(f,2*msignalsdft);
    zoom xon;
    zoom(300);
    grid on
    title('DFT of Stimulus');
    xlabel('Frequency (Hz)'); 
    ylabel('Amplitude (dB)');
    subplot(2,2,2)
    plot(((0:(n-1)) * dt),signal_r)
    title('Response');
    xlabel('Time (s)'); 
    ylabel('Amplitude (degrees)');
    subplot(2,2,4)  
    plot(f,2*msignalrdft);
    zoom xon;
    zoom(300);
    grid on
    title('DFT of Response');
    xlabel('Frequency (Hz)'); 
    ylabel('Amplitude (dB)');

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 6];
    img_name = strcat(strcat(folder,'\'), strcat('DFT_',name));
    print ('-dpng', img_name);

    %{
    pause('on');
    pause(1);
    %}
    DFT_s = [f', 2*msignalsdft];
    DFT_r = [f', 2*msignalrdft];    
    csvwrite(strcat(strcat(folder,'\'),strcat(name, '_dft.csv')), DFT_s);
    csvwrite(strcat(strcat(folder,'\'),strcat(name, '_dft.csv')), DFT_r);
    close(gcf);

end

function [signaldft_HannWnd] = Signal_Mod(signal, L, nfft, padding)
    signal_HannWnd = signal.*hanning(L, 'periodic'); %Apply Hann Window
    signal_HannWnd = padarray(signal_HannWnd, padding); %zero padding of signal
    signaldft_HannWnd = abs(fftshift(fft(signal_HannWnd,nfft))); %Absolute values of FFT of signal
    return
end