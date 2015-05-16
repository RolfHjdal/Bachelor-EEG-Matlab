close all; clear all;


%last alle EEGData
%load EEGData.mat
% 
% load neutral.mat
% neutral_1 = neutral.signals.values(1:1280, from:to);
% load neutral_2.mat
% neutral_2 = neutral_2.signals.values(1:1280, from:to);
% load neutral_3.mat
% neutral_3 = neutral_3.signals.values(1:1280, from:to);
% load push.mat
% push_1 = ScopeData.signals.values(1:1280, from:to);
% load push_2.mat
% push_2 = ScopeData2.signals.values(1:1280, from:to);
% load push_3.mat

S= load('newAssburgers.mat')
for f = fieldnames(S)'
   disp(['Field named: ' f{1} ]);
   disp('Has value ')
   disp(S.(f{1}));
   
end

return
from=2;
to=15;
N=1024;
fs=128;
input = zeros(7,200);
target = [ones(1,100) zeros(1,100) ; zeros(1,100) ones(1,100)]

%%For kvar matrise; finn differanse mellom relevante elektroder
for i = 1:7
    neutral_1_diff(:,i) = neutral_1(:,i)-neutral_1(:,15-i);
end

%remove DC from all channels
for i=1:7
neutral_1_diff(:,i)=neutral_1_diff(:,i)-mean(neutral_1_diff(:,i));
end

%FFT
% Frekvensoppløsninga delta(f) = (1/N)*fs = 128/1024 = 0,125Hz
FFT_neutral_1 = abs(fft(neutral_1_diff, N));
stem(FFT_neutral_1(1:N/2, 1));
ylabel('|x(n)|')
xlabel('Sample n')
title('delta(f) = fs/N = 1Hz')
