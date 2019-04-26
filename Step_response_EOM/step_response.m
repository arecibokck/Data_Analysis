function plotScopeTrace
    filename = 'E:\Master Thesis\Data\Step_Response_EOM\02-25-19\maxpgain2Mhzsigain500KHzdiffgain\scope_trace_average.csv';
    delimiter = ',';
    startRow = 3;
    formatSpec = '%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    scopetrace = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);

    time = scopetrace{:, 1} + 5.00E-05; ...- 0.01098406;
    step = smooth(normalize(scopetrace{:, 2}));
    in_loop = normalize(scopetrace{:, 3});
    error = scopetrace{:, 4};
    
    csvwrite('step_data.csv', [time step]);
    csvwrite('response_data.csv', [time in_loop]);
    
    hold on;
    plot(time, step);
    plot(time, in_loop);
    plot(time, error);
    hold off;
    
    axis tight;
    %title('Step Response of Analog Laser Servo Feedback Loop with EOM');
    title('Averaged Step Response');
    legend({'Step', 'In-Loop PD', 'Error'}, 'Location', 'west');
    xlabel('Time(s)') 
    ylabel('Normalized Amplitude')
   
end

function normalized_vals = normalize(vals)
    normalized_vals = (vals - min(vals)) / (max(vals) - min(vals)); 
end
