close all; clear all;

input = zeros(7,200); %Preallocate
target = [ones(1,100) zeros(1,100) ; zeros(1,100) ones(1,100)]; %Build target vector
from=3;
to=16;
N=1024;
fs=128;
inputFrequency=15; %Frequency to pass to the neural network
inputFrequency2=20; %Frequency to pass to the neural network
inputFrequency3=25; %Frequency to pass to the neural network
inputFrequency4=30; %Frequency to pass to the neural network
inputFrequency5=18; %Frequency to pass to the neural network
inputFrequency6=13; %Frequency to pass to the neural network
k=inputFrequency*(N/fs); %Find corresponding point k in fft to frequency
k2=inputFrequency2*(N/fs); %Find corresponding point k in fft to frequency
k3=inputFrequency3*(N/fs); %Find corresponding point k in fft to frequency
k4=inputFrequency4*(N/fs); %Find corresponding point k in fft to frequency
k5=inputFrequency5*(N/fs); %Find corresponding point k in fft to frequency
k6=inputFrequency6*(N/fs); %Find corresponding point k in fft to frequency
deltaL = 1 * (N/fs);
deltaH = 4 * (N/fs);
thetaL = 4 * (N/fs);
thetaH = 7 * (N/fs);
alphaL = 7 * (N/fs);
alphaH = 15 * (N/fs);
betaL1 = 15 * (N/fs);
betaH1 = 20 * (N/fs);
betaL2 = 20 * (N/fs);
betaH2 = 25 * (N/fs);
betaL3 = 25 * (N/fs);
betaH3 = 30 * (N/fs);
currentDataSet=1; %Dataset currently being worked on
S= load('matlab.mat');

%Loop through all elements
for f = fieldnames(S)'
    channelRAW=S.(f{1})(1:1280, from:to);

%%Find potential between electrode pairs
for i = 1:7
    channel(:,i) = channelRAW(:,i)-channelRAW(:,15-i); 
end

%Remove DC from all channels by subtracting mean
for i=1:7
    channel(:,i)=channel(:,i)-mean(channel(:,i));
end

% FFT
FFTchannel = abs(fft(channel, N));

%Build input vector for the neural network
for i=1:7
    input(i, currentDataSet) = FFTchannel(k, i);
    inputMean(i, currentDataSet) = mean(FFTchannel(deltaL:deltaH, i));
end

%Build input vector for the neural network
for i=8:1:14
    input(i, currentDataSet) = FFTchannel(k2, i-7);
    inputMean(i, currentDataSet) = mean(FFTchannel(thetaL:thetaH, i-7));
end

%Build input vector for the neural network
for i=15:1:21
    input(i, currentDataSet) = FFTchannel(k3, i-14);
    inputMean(i, currentDataSet) = mean(FFTchannel(alphaL:alphaH, i-14));
end

%Build input vector for the neural network
for i=22:1:28
    input(i, currentDataSet) = FFTchannel(k4, i-21);
    inputMean(i, currentDataSet) = mean(FFTchannel(betaL1:betaH1, i-21));
end

%Build input vector for the neural network
for i=29:1:35
    input(i, currentDataSet) = FFTchannel(k5, i-28);
    inputMean(i, currentDataSet) = mean(FFTchannel(betaL2:betaH2, i-28));
end

%Build input vector for the neural network
for i=36:1:42
    input(i, currentDataSet) = FFTchannel(k6, i-35);
    inputMean(i, currentDataSet) = mean(FFTchannel(betaL3:betaH3, i-35));
end
currentDataSet=currentDataSet+1;
end
%% Plot illustration of the final ANN input matrix 42x200 with mean values
mean1=mean(FFTchannel(deltaL:deltaH, 1));
mean2=mean(FFTchannel(thetaL:thetaH, 1));
mean3=mean(FFTchannel(alphaL:alphaH, 1));
mean4=mean(FFTchannel(betaL1:betaH1, 1));
mean5=mean(FFTchannel(betaL2:betaH2, 1));
mean6=mean(FFTchannel(betaL3:betaH3, 1));
allMean=[mean1 mean2 mean3 mean4 mean5 mean6];

k=inputFrequency*(N/fs); %Find corresponding point k in fft to frequency
k2=inputFrequency2*(N/fs); %Find corresponding point k in fft to frequency
k3=inputFrequency3*(N/fs); %Find corresponding point k in fft to frequency
k4=inputFrequency4*(N/fs); %Find corresponding point k in fft to frequency
k5=inputFrequency5*(N/fs); %Find corresponding point k in fft to frequency
k6=inputFrequency6*(N/fs); %Find corresponding point k in fft to frequency
close all; 
%plots the FFT together with mean of frequency bands.
%repmat(mean2, thetaH-thetaL+1)
area(1:512, FFTchannel(1:512,1), 'LineWidth', 0.5); hold on
plot(1:deltaH, ones(1, deltaH)*mean1, 'g', 'LineWidth', 2); 
plot(thetaL:thetaH, ones(1, thetaH-thetaL+1)*mean2, 'c', 'LineWidth', 2)
plot(alphaL:alphaH,  ones(1, alphaH-alphaL+1)*mean3, 'm', 'LineWidth', 2)
plot(betaL1:betaH1,  ones(1, betaH1-betaL1+1)*mean4, 'r', 'LineWidth', 2)
plot(betaL2:betaH2,  ones(1, betaH2-betaL2+1)*mean5,  'r',   'LineWidth', 2)
plot(betaL3:betaH3,  ones(1, betaH3-betaL3+1)*mean6, 'r', 'LineWidth', 2); hold off
title('FFT and mean values (Channel 1)')
ylabel('Magnitude |X(k)|')
xlabel('Datapoint k in FFT')
h=legend('FFT', '\delta (1-4Hz)', '\theta (4-7Hz)', '\alpha (7-15Hz)', '\beta Low (15-20Hz)', '\beta Medium (20-25Hz)', '\beta High (25-30Hz)')

grid on
axis([0 512 0 3500])
%% illustration of 7x200 input matrix @20Hz
mean1=mean(FFTchannel(deltaL:deltaH, 1));
mean2=mean(FFTchannel(thetaL:thetaH, 1));
mean3=mean(FFTchannel(alphaL:alphaH, 1));
mean4=mean(FFTchannel(betaL1:betaH1, 1));
mean5=mean(FFTchannel(betaL2:betaH2, 1));
mean6=mean(FFTchannel(betaL3:betaH3, 1));
allMean=[mean1 mean2 mean3 mean4 mean5 mean6];
k=inputFrequency*(N/fs); %Find corresponding point k in fft to frequency
k2=inputFrequency2*(N/fs); %Find corresponding point k in fft to frequency
k3=inputFrequency3*(N/fs); %Find corresponding point k in fft to frequency
k4=inputFrequency4*(N/fs); %Find corresponding point k in fft to frequency
k5=inputFrequency5*(N/fs); %Find corresponding point k in fft to frequency
k6=inputFrequency6*(N/fs); %Find corresponding point k in fft to frequency
close all; 
%plots the FFT together with mean of frequency bands.
%repmat(mean2, thetaH-thetaL+1)
area(1:512, FFTchannel(1:512,1), 'LineWidth', 0.5); hold on
%plot(1:deltaH, ones(1, deltaH)*mean1, 'g', 'LineWidth', 2); 
%plot(thetaL:thetaH, ones(1, thetaH-thetaL+1)*mean2, 'c', 'LineWidth', 2)
%plot(alphaL:alphaH,  ones(1, alphaH-alphaL+1)*mean3, 'm', 'LineWidth', 2)
%plot(betaL1:betaH1,  ones(1, betaH1-betaL1+1)*mean4, 'r', 'LineWidth', 2)
%plot(betaL2:betaH2,  ones(1, betaH2-betaL2+1)*mean5,  'r',   'LineWidth', 2)
plot(betaL3:betaL3,  ones(1, betaH3-betaL3+1)*mean6, 'r', 'LineWidth', 2); hold off
title('FFT (Channel 1)')
ylabel('Magnitude |X(k)|')
xlabel('Datapoint k in FFT')
h=legend('FFT', '\beta 20Hz')
line([k2 k2], [0.1 750], 'Color', 'r', 'LineWidth', 2) %plot marker @20Hz
grid on
axis([0 512 0 3500])
%% illustration of 7x200 input matrix using several frequencies
mean1=mean(FFTchannel(deltaL:deltaH, 1));
mean2=mean(FFTchannel(thetaL:thetaH, 1));
mean3=mean(FFTchannel(alphaL:alphaH, 1));
mean4=mean(FFTchannel(betaL1:betaH1, 1));
mean5=mean(FFTchannel(betaL2:betaH2, 1));
mean6=mean(FFTchannel(betaL3:betaH3, 1));
allMean=[mean1 mean2 mean3 mean4 mean5 mean6];
k=8*(N/fs); %Find corresponding point k in fft to frequency
k2=13*(N/fs); %Find corresponding point k in fft to frequency
k3=16*(N/fs); %Find corresponding point k in fft to frequency
k4=20*(N/fs); %Find corresponding point k in fft to frequency
k5=25*(N/fs); %Find corresponding point k in fft to frequency
k6=30*(N/fs); %Find corresponding point k in fft to frequency
close all; 
%plots the FFT together with mean of frequency bands.
%repmat(mean2, thetaH-thetaL+1)
area(1:512, FFTchannel(1:512,1), 'LineWidth', 0.5); hold on
plot(1:1,  ones(1, 1)*mean1, 'green', 'LineWidth',  2) 
plot(thetaL:thetaL,  ones(1, thetaH-thetaH+1)*mean2, 'cyan', 'LineWidth', 2)
plot(alphaL:alphaL,  ones(1, alphaH-alphaH+1)*mean3, 'magenta', 'LineWidth',  2)
plot(betaL1:betaL1,  ones(1, betaH1-betaH1+1)*mean4, 'red', 'LineWidth',  2)
plot(betaL2:betaL2,  ones(1, betaH2-betaH2+1)*mean5, 'red', 'LineWidth',  2)
plot(betaL3:betaL3,  ones(1, betaH3-betaH3+1)*mean6, 'r', 'LineWidth',  2); hold off
title('FFT (Channel 1)')
ylabel('Magnitude |X(k)|')
xlabel('Datapoint k in FFT')
legend('FFT', '\alpha 8Hz', '\alpha 13Hz', '\beta 16Hz', '\beta 20Hz', '\beta 25Hz', '\beta 30Hz')

line([k k], [0.1 750], 'Color', 'g', 'LineWidth', 2) %plot marker @8Hz
line([k2 k2], [0.1 750], 'Color', 'c', 'LineWidth', 2) %plot marker @13Hz
line([k3 k3], [0.1 750], 'Color', 'm', 'LineWidth', 2) %plot marker @16Hz
line([k4 k4], [0.1 750], 'Color', 'r', 'LineWidth', 2) %plot marker @20Hz
line([k5 k5], [0.1 750], 'Color', 'r', 'LineWidth', 2) %plot marker @25Hz
line([k6 k6], [0.1 750], 'Color', 'r', 'LineWidth', 2) %plot marker @30Hz
grid on
axis([0 512 0 3500])
%%
 %Plot data
 %mean values for single plot
for j=1:200
inputPlot(j)=mean(input(:,j));
end

%%

figure
%A=4, B=7
channelA=4;
channelB=7;

scatter(input(channelA,1:100), input(channelB,1:100), 'o'); hold on %meditation
scatter(input(channelA, 101:200), input(channelB,101:200), '+'); %push
hold off
grid on
title('Scatter Push/Meditation')
xlabel('Channel A')
ylabel('Channel B')
legend('Push','Meditation')
%saveas(figure(1), 'Scatter 6-72.png')
axis([0 1000 0 1000])
%%
figure
channelA=6;
channelB=7;

scatter(input(channelA,1:100), input(channelB,1:100), 'o'); hold on %meditation
scatter(input(channelA, 101:200), input(channelB,101:200), '+'); %push
hold off
grid on
title('Scatter Push/Meditation')
xlabel('Channel A')
ylabel('Channel B')
legend('Push','Meditation')
axis([0 1000 0 1000])

figure
channelA=1;
channelB=3;

scatter(input(channelA,1:100), input(channelB,1:100), 'o'); hold on %meditation
scatter(input(channelA, 101:200), input(channelB,101:200), '+'); %push
hold off
grid on
title('Scatter Push/Meditation')
xlabel('Channel A')
ylabel('Channel B')
legend('Push','Meditation')
axis([0 1000 0 1000])