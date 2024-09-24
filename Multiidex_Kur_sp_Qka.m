%% multi-index 
% ks=kurtosis(y)/skewness(y); %峭度偏斜度比
% peakingfactor(y)%峰值因子
% sparsity(y)%稀疏度

%% ini
clear
clc
close all
%% Two fault modes with common resonant frequency
fr1=3;  % convolution freq.
a1=1;    % amplitude for sig1;
fn1=625; % resonant freq.
zeta=400; % decaying rate
T1=1/16;  % fault frequency for sig1;
fs=10000;  % sampling freq.

% Generating the signal
t1=0:1/fs:T1;
sig_1=(a1*exp(-zeta*t1).*cos(2*pi*fn1*t1));  % single impulse

t=0:1/fs:2;
sig0=zeros(size(t));
for i=1:fix(t(end)/T1)
    if i==1 || i==fix(t(end)/T1)
        inde=max(find(t<=T1*(i-1)));
    else
        inde=max(find(t<=T1*(i-1)+(rand(1)-0.5)*2*1/100*T1));
%         inde=max(find(t<=T1*(i-1)+dd(mod(round(rand(1)*10),3)+1)*13/100*T1));
%         inde=max(find(t<=T1*(i-1)+dd(randperm(3,1))*0.75/100*T1));
    end
    sig0(inde:inde+length(sig_1)-1)=sig_1;
end
% sig1=(cos(2*pi*fr1*t)+1).*sig0;  % fault signal 1
sig1=sig0;
sig=awgn(sig1,-5,'measured');
% figure,plot(t,sig,'b')
figure,plot(t,sig,'color',0.5*[1 1 1])
hold on
plot(t,sig1,'b-')
axis tight
xlim([0 2])
xlabel('time [s]'),ylabel('Amplitude')
% frequency analysis
pp=abs(fft(sig))/length(sig)*2;
pp(1)=0;
ff=(0:length(pp)-1)/length(pp)*fs;
n1=floor(length(pp)/2);
%figure,plot(ff(1:n1),pp(1:n1),'b')
% xlim([400 1000])
xlabel('Frequency [Hz]'),ylabel('Magnitude')
y=sig;
y=y(1:20000);
soft = @(x, T) max(1 - T./abs(x), 0) .* x;

%% signal parameters
fs=10000;  % sampling freq.
t=(0:length(y)-1)/fs;
figure,plot(t,y,'b')





%% multi-index matrix TQWT  (Q,J,ka inter)
gamma=0.85;
Q_range = 1:1:10;               
r=3;
ka_range= 0.3:0.1:1;
N=200;
matrix=cell(length(Q_range),length(ka_range));
f=1/T1;
delt_p=1.5;
for i =1:length(Q_range)
    Q=Q_range(i);
    J=min(10,computeJmax(N,Q,r));
    for j=1:length(ka_range)
        ka=ka_range(j);
        [w,~]=TQWT_SR_GMC_penalty_fun(y,Q,r,J,gamma,ka,0);
        y_GMC = itqwt(w,Q,r,length(y));
        matrix{i,j}=[kurtSkew(y_GMC);ES_EHNR(y_GMC);sparsity(y_GMC)]; %获得了参数矩阵
    end
end

kurtSkew_vector=zeros(0);
ES_EHNR_vector=zeros(0);
sparsity_vector=zeros(0);
 for i=1:length(Q_range)
    for j=1:length(ka_range)
        kurtSkew_vector(end+1)=matrix{i,j}(1);
        ES_EHNR_vector(end+1)=matrix{i,j}(2);
        sparsity_vector(end+1)=matrix{i,j}(3);
    end
 end
matrix_index=[kurtSkew_vector;ES_EHNR_vector;sparsity_vector];


%% choose Q,ka
[Fj,loc] = chooseband(matrix_index);
% 现在要找到1oc在10×7这个矩阵中的位置，并返回下标 即Q的值和ka的值 注意Q和ka的值最后要换算到Q、ka_range中
Q_loc=ceil(loc/length(ka_range));%参数矩阵横坐标 
ka_loc=loc-(Q_loc-1)*length(ka_range);%参数矩阵纵坐标
Q=Q_range(Q_loc);
ka=ka_range(ka_loc);
J=min(10,computeJmax(N,Q,r));
%% TQWT sparse representation
[x,v]=TQWT_SR_GMC_penalty_fun(y,Q,r,J,gamma,ka,0);
y_GMC = itqwt(x,Q,r,length(y));
figure,plot(t,y_GMC)
xlim([0,0.4])
figure,plot(t,sig1(1:length(y)),t,y_GMC,'r--')
legend('故障特征信号','去噪信号')
xlabel('时间[s]'),ylabel('幅值')

figure,plot(t,sig1(1:length(y)),t,y_GMC,'r--')
xlim([0.1,0.15])
legend('故障特征成分','小波基函数')
xlabel('时间[s]'),ylabel('幅值')
set(gca,'xticklabel',[]),set(gca,'yticklabel',[])
% figure,plot(t,y_GMC)
% xlabel('时间[s]'),ylabel('幅值')
% xlim([0,1])
xlim([4000,6000])
xlabel('采样时间')
ylabel('幅值')
set(gca,'xticklabel',[]),set(gca,'yticklabel',[])
title('Q=1,r=3')


figure,plot(sig_1)
xlabel('时间[s]')
ylabel('幅值')
xlim([0,300])


blp=abs(fft(abs(hilbert(y_GMC))))/length(y_GMC)*2;
figure,plot(blp)

blp(1)=0;
pl=(0:length(y_GMC)-1)/length(y_GMC)*fs;
figure,plot(pl(1:round(length(y_GMC)/2)),blp(1:round(length(y_GMC)/2)))
xlabel('频率 [Hz]'),ylabel('幅值')
xlim([0,350])