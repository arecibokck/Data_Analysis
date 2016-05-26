function Correlation(s, r, Fs, folder, name)
    
    N = length(s);
    
    signal_s = detrend(s); %detrend data to remove DC Mean Component (subtract average from each 
    signal_r = detrend(r); %data-point) - Now there is no peak at 0 Hz 

    [y_a, lag_a] = xcorr(signal_s, 'coeff')
    [~,A] = max(abs(y_a));
    lagDiff_a = lag_a(A)
    timeDiff_a = lagDiff_a/Fs
    y_a = y_a(N/1.5:N-1, 1);
    t_a = (-length(y_a)/2:length(y_a)/2-1)/Fs


    [y_c, lag_c] = xcorr(signal_s, signal_r, 'coeff')
    [~,C] = max(abs(y_c));
    lagDiff_c = lag_c(C)
    timeDiff_c = lagDiff_c/Fs
    y_c = y_c(N/1.5:N-1,1);
    t_c = (-length(y_c)/2:length(y_c)/2-1)/Fs
    figure(2)
    set(0,'defaultaxesFontName', 'CMU Serif Roman')
    set(0,'defaultaxesFontSize', 12)
    plot(t_a, y_a, 'b.');
    hold on;
    plot(t_c, y_c, 'r.');
    title('Correlation');
    xlabel('Lag (s)'); 
    ylabel('Normalized CF');
    legend('AutoCorr \theta_{Body-Body}', 'CrossCorr \theta_{Head-Body}', 'Location', 'SouthEast');
    hold off;
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 6];
    img_name = strcat(folder, 'Correlation_', name);
    print ('-dpng', img_name);
    %{
    pause('on');
    pause(1);
    %}
    close(gcf);
end