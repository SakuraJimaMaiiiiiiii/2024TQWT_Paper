%% Slip_Simulated_Signal_Analysis

clear
close all
clc

%% chosen signal
% original signal
M = csvread('468.csv',1,0);
sig0=M(:,1);
figure,plot(sig0,'y--')
fs=25600;  % sampling freq.

% 20000length signal
y = sig0;
y = y(1:20000);
figure,plot(y,'r--')
xlabel('time [s]'),ylabel('Amplitude')
axis tight
blp = abs(y);
figure,plot(blp)
t=(0:length(y)-1)/fs;
figure,plot(t,y,'b')
xlabel('time [s]'),ylabel('Amplitude')

% guzhang signal
t_end= 1.28;
tt = length(sig0)/t_end/fs; 
guzhang_signal = sig0(1:tt:length(sig0));
t=(0:length(y)-1)/fs;
figure,plot(t,guzhang_signal(1:length(t)),'b--')
x=guzhang_signal;

%% TQWT denoising
x=x(1:20000);
xzero=zeros(size(x));
N=length(x);
Q = 3; r = 3; J = 15; % High Q-factor wavelet transform
Jtext=computeJmax(N,Q,r);
J=min(15,Jtext);
[wlets, now] = ComputeWavelets(N,Q,r,J,'radix2');
w = tqwt(x,Q,r,J);              % TQWT
wzero = tqwt(xzero,Q,r,J);      %follow sig

norm_w=w{1}/now(1);      %guiyihua
for i=2:length(now)
    norm_w=[norm_w(:);w{i}(:)/now(i)];
end

%% definition of the threshold
thres=mean(norm_w)+std(norm_w)*2.5;
% thres=thselect(norm_w,'minimaxi'); % 'heursure','sqtwolog','minimaxi','rigrsure'
%% wavelet coefficient 
alpha=0.9;
for i=1:J+1
    wzero{i}=compute_soft(w{i},thres,alpha);
end
%% inverse TQWT and obtain the denoised result
y = itqwt(wzero,Q,r,length(x));  
                                  %figure,plot(y)      reconstrust sig
                                  %figure,plot(sig0(1:length(y)),'b')
figure,plot(x,'b-')
hold on 
plot(y,'r')
xlabel('times [s]')
ylabel('幅值')
legend('468时刻故障信号','去噪信号')


f = linspace(0,fs/2,N/2);
A1 = abs(y)/(N/2);
figure,plot(f,A1(1:N/2))


blp=abs(fft(abs(hilbert(y))))/length(y)*2;
blp(1)=0;
pl=(0:length(y)-1)/length(y)*fs;
figure,plot(pl(1:round(length(y)/2)),blp(1:round(length(y)/2)),'r')
xlabel('frequency [HZ]')
ylabel('Magnitude')

%外圈
% line([115.6 115.6],[0 0.8]);
% hold on
% line([231.2 231.2],[0 0.8]);
% hold on
% line([346.8 346.8],[0 0.8]);
% hold on
% line([462.4 462.4],[0 0.8]);
% hold on
% line([578 578],[0 0.8]);

%内圈
% line([184.4 184.4 ],[0 0.8]);
% hold on
% line([368.8 368.8],[0 0.8]);
% hold on
% line([553.2 553.2],[0 0.8]);
% hold on
% line([737.6 737.6],[0 0.8]);
% hold on
% line([922 922],[0 0.8]);
% 
