function data = pre_process_data_v1(dat1)

k = 250; % point of transition from linear to log
data = dat1;
data(data>k) = k*(log(data(data>k))-(log(k)-1));
data(data<-k) = -k*(log(-data(data<-k))-(log(k)-1));

% notch filter 50Hz
fc = 50; fs = 256;
load('50Hz_notch.mat'); % Num and Den out here
%[num,den] = iirnotch(2*fc/fs,4/fs);
data = filter(Num, Den, data')';

% High pass filter
[B,A] = butter(4, 1/256, 'high'); % additional attenuation at 50Hz, 3dB~32Hz
data = filter(B, A, data')';
