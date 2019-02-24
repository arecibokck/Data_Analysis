% script to plot Bode plots.
 
% step=importdata('data_500_1.csv');                       %data for step function
% response= importdata('data_500_2.csv');                  %data for response function
 
step = importdata('step_data.csv'); 
response = importdata('response_data.csv');   ... importdata('response_data.csv');                  

% step=step(80:end,:);                                     %delete the first 80 entries 
% response=response(80:end,:);

periodStep = (step(end, 1) - step(1, 1));                %sampling period
samplingTimeStep = (step(2, 1) - step(1, 1));            %time between two data points

%derivativeStep =[.5*(step(2:end, 1) + step(1:(end-1),1)), ((step(2:end,2) - step(1:(end-1),2))/(samplingTimeStep))];    

periodResponse = (response(end, 1) - response(1, 1));
samplingTimeResponse = (response(2, 1) - response(1, 1));

derivativeResponse =[.5*(response(2:end, 1) + response(1:(end-1),1)), ((response(2:end,2) - response(1:(end-1),2))/(samplingTimeResponse))];  %derivative of the response curve (impulse function)
derivativeStep = [.5*(step(2:end, 1) + step(1:(end-1),1)), ((step(2:end,2) - step(1:(end-1),2))/(samplingTimeStep))];

response(:,2)=response(:,2)-mean(response(1:100,2));      % make offset zero
norm= mean(response((end-200):end,2));                    % choose last 200 data points to normalize the curve
response(:,2)=response(:,2)/norm;                         %normalize to unit height

trigger=1;
while step(trigger+1,2)-step(trigger,2)<=0.0009 %find the trigger
    trigger=trigger+1;
end

response(:,1)=response(:,1)-response(trigger,1);
step(:,1)=step(:,1)-step(trigger,1);                      %trigger begins at t=0

%paddedResponse= cat(1,derivativeResponse,paddingData);

% responseDerivative=cat(1,derivativeResponse,[
subplot(2,2,1);
hold on
plot(response(:,1).*10^6,response(:,2));
plot(step(:,1).*10^6,step(:,2),'r');
xlabel('time(\mus)');
ylabel('step response function');
grid on
hold off
axis tight;


subplot(2,2,2);
hold on
plot(derivativeResponse(:,1).*10^6,derivativeResponse(:,2).*10^-6,'g');
% plot(derivativeStep(:,1).*10^6,derivativeStep(:,2).*10^-6,'g');
% plot(derivativeResponse(:,1).*10^6,smooth(derivativeResponse(:,2).*10^-6,3))
xlabel('time(\mus)');
ylabel('Transfer function (1/ \mus)');
hold off
axis tight;




fs = 1/samplingTimeResponse;                                         % Frequency span
resDeriv=derivativeResponse(trigger:end,2);                          % take values from the trigger onwards
m = length(resDeriv);                                                % Window length
n = m*41;                                                            %zero padding
fourierResponse=fft(resDeriv,n);
p=fourierResponse((n/2):end);
c=fourierResponse(2:(n/2-1));
y = cat(1,p,c);

% f = (-((n/2)-1):((n/2)-1))*(fs/m);                                    % Frequency range from negative to positive


freqRange= (0:((n/2)-3))*(fs/n);
% fa = (-((n/2)):((n/2)-1))*(fs/n);
phase=unwrap(angle(c));
spectralPower= (abs(c)/fs).^2;
bode=10*log10(spectralPower);


subplot(2,2,3);

% plot(f,((abs(y)).^2))
semilogx(freqRange,bode);
k=1;
while bode(k)>-20
    k=k+1;
end

xlim([0 freqRange(k)]);
b=1;
while bode(b)>-3
    b=b+1;
end
xlabel('Frequency (Hz)')
ylabel('Spectral Power(dB)')
xl = get(gca, 'xlim')
text(900,-18,cat(2,'Bandwidth: ',[num2str(round(freqRange(b)/1000)),' kHz']));



subplot(2,2,4);
phaseDegree=(phase/pi)*180;
semilogx(freqRange,phaseDegree,'g')
q=1;
while phaseDegree(q)>-180
    q=q+1;
end
ylim([-250 0]);
xlabel('Frequency (Hz)')
ylabel('Phase(degrees)')
text(150,-225,cat(2,'Cut-off frequency: ',[num2str(round(freqRange(q)/1000)),' kHz']));


hold off



