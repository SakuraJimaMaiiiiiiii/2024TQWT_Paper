%% multi-index 
% ks=kurtosis(y)/skewness(y); %峭度偏斜度比
% ES_EHNR(y)%谐波与噪声能量比
% sparsity(y)%稀疏度


%% import signal
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
% figure,plot(t,sig1,'color',0.5*[1 1 1])
% xlim([0,1]);
% xlabel('时间[s]'),ylabel('幅值A / V')
sig=awgn(sig1,-20,'measured');
%figure,plot(t,sig,'b')
figure,plot(t,sig,'color',0.5*[1 1 1])
hold on
plot(t,sig1,'b-')
axis tight
xlim([0 1])
xlabel('time [s]'),ylabel('Amplitude')
% frequency analysis
pp=abs(fft(sig))/length(sig)*2;
pp(1)=0;
ff=(0:length(pp)-1)/length(pp)*fs;
n1=floor(length(pp)/2);
figure,plot(ff(1:n1),pp(1:n1),'color',0.5*[1 1 1])
% xlim([0 800])
xlabel('频率 [Hz]'),ylabel('幅值e')
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
f=1/T1;
delt_p=1.5;
%确定J的范围 用以初始化元胞数组
maxJ=zeros();
for i=1:length(Q_range)
    maxJ(i)=computeJmax(N,i,r);
end
MAXJ=max(maxJ);

matrix=cell(length(Q_range),MAXJ,length(ka_range));

for i =1:length(Q_range)
        Q=Q_range(i);
        J_range=computeJmax(N,Q,r);
        J_len=1:J_range;
    for j=1:MAXJ
            for k=1:length(ka_range)
                if j<=J_range
                    J=J_len(j);
                    ka=ka_range(k);
                    [w,~]=TQWT_SR_GMC_penalty_fun(y,Q,r,J,gamma,ka,0);
                    y_GMC = itqwt(w,Q,r,length(y));
                    matrix{i,j,k}=[kurtSkew(y_GMC);ES_EHNR(y_GMC);sparsity(y_GMC)]; %获得了参数矩阵
                else
                    matrix{i,j,k}=[0;0;0];
                end
            end
    end
end

kurtSkew_vector=zeros(0);
ES_EHNR_vector=zeros(0);
sparsity_vector=zeros(0);
 for i=1:size(matrix,1)
    for j=1:size(matrix,2) 
        for k=1:size(matrix,3)
            kurtSkew_vector(end+1)=matrix{i,j,k}(1);
            ES_EHNR_vector(end+1)=matrix{i,j,k}(2);
            sparsity_vector(end+1)=matrix{i,j,k}(3);
        end
    end
 end  
matrix_index=[kurtSkew_vector;ES_EHNR_vector;sparsity_vector];


%% choose Q,ka
[Fj,loc] = chooseband(matrix_index);
% 现在要找到1oc在三维矩阵中的位置，由于choosebang的结果是一行向量，即
%由i,j,k按读取顺序排列的 所以需要进行操作访问下标 注意Q和ka的值最后要换算到Q、ka_range中，J则还与Q有关
Q_loc=ceil(loc/size(matrix,2)/size(matrix,3));%参数矩阵第一层坐标 Q
J_loc=ceil((loc-size(matrix,2)*size(matrix,3)*(Q_loc-1))/size(matrix,3));%参数矩阵第二层坐标 Q
%J的最终值还需要找到Q_loc决定的J_range
ka_loc=loc-(Q_loc-1)*size(matrix,2)*size(matrix,3)-(J_loc-1)*size(matrix,3);%参数第三层坐标 ka

Q_choose=Q_range(Q_loc);
%J的最终值还需要找到Q_loc决定的J_ranged 然后再找J的位置
J_vetor=1:computeJmax(N,Q_loc,r);
J_fin=J_vetor(J_loc);
ka_choose=ka_range(ka_loc);

%% parm sort
d=imag(Fj);
figure,plot(sort(d,'descend')) ;
ylabel('融合特征指标')
xlabel('数据长度')


%% TQWT sparse representation
[x,v]=TQWT_SR_GMC_penalty_fun(y,Q_choose,r,J_fin,gamma,ka_choose,0);
y_GMC = itqwt(x,Q_choose,r,length(y));
figure,plot(t,y_GMC,'color',0.5*[1 1 1])
plot(t,sig1(1:length(y)))
xlim([0,0.4])
figure,plot(t,sig1(1:length(y)),t,y_GMC,'r--')
legend('故障特征信号','去噪信号')
xlabel('时间[s]'),ylabel('幅值')

figure,plot(t,sig1(1:length(y)),t,y_GMC,'r--')
xlim([0.12,0.135])
legend('故障特征成分','小波基函数')
xlabel('时间'),ylabel('幅值A/ V')
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
xlim([0,250])


blp=abs(fft(abs(hilbert(y_GMC))))/length(y_GMC)*2;
figure,plot(blp,'color',0.5*[1 1 1])

blp(1)=0;
pl=(0:length(y_GMC)-1)/length(y_GMC)*fs;
figure,plot(pl(1:round(length(y_GMC)/2)),blp(1:round(length(y_GMC)/2)),'color',0.5*[1 1 1])
xlabel('频率 [Hz]'),ylabel('幅值')
xlim([0,350])




