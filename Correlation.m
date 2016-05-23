function Correlation(stimulus, response, Fs, folder, name)
    
    signal_s = detrend(csvread(stimulus));
    signal_r = detrend(csvread(response));

    [y_a, lag_a] = xcorr(signal_s, 'coeff')
    [~,A] = max(abs(y_a));
    lagDiff_a = lag_a(A)
    timeDiff_a = lagDiff_a/Fs
    y_a = y_a(10000:11000,1);
    t_a = (-length(y_a)/2:length(y_a)/2-1)/Fs


    [y_c, lag_c] = xcorr(signal_s, signal_r, 'coeff')
    [~,C] = max(abs(y_c));
    lagDiff_c = lag_c(C)
    timeDiff_c = lagDiff_c/Fs
    y_c = y_c(10000:11000,1);
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
    img_name = strcat(strcat(folder,'\'), strcat('Correlation_',name));
    print ('-dpng', img_name);
    %{
    pause('on');
    pause(1);
    %}
    close(gcf);
end