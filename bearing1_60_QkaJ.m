%% Slip_Simul ated_Signal_Analysis

clear
close all
clc

%% chosen signal
% original signal
M = csvread('60.csv',1,0);
sig0=M(:,1);
figure,plot(sig0,'y--')
fs=26500;  % sampling freq.

% 20000length signal
y = sig0;
y = y(1:20000);
figure,plot(y)
xlabel('时间 [s]'),ylabel('幅值')
xlim([1*10^4,1.5*10^4])
xticklabels({'1','1.05','1.1','1.15','1.2','1.25','1.3','1.35','1.4','1.45','1.5'})


%% multi-index matrix TQWT  (Q,J,ka inter)
gamma=0.85;
Q_range = 1:1:10;               
r=3;
ka_range= 0.3:0.1:1;
N=200;

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
                    [w,~]=TQWT_SR_GMC_penalty_fun(y',Q,r,J,gamma,ka,0);
                    y_GMC = itqwt(w,Q,r,length(y));
                    matrix{i,j,k}=[kurtosis(y_GMC);peakingfactor(y_GMC);sparsity(y_GMC)]; %获得了参数矩阵
                else
                    matrix{i,j,k}=[0;0;0];
                end
            end
    end
end

kurtosis_vector=zeros(0);
peakingfactor_vector=zeros(0);
sparsity_vector=zeros(0);
 for i=1:size(matrix,1)
    for j=1:size(matrix,2)
        for k=1:size(matrix,3)
            kurtosis_vector(end+1)=matrix{i,j,k}(1);
            peakingfactor_vector(end+1)=matrix{i,j,k}(2);
            sparsity_vector(end+1)=matrix{i,j,k}(3);
        end
    end
 end
matrix_index=[kurtosis_vector;peakingfactor_vector;sparsity_vector];


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



%% TQWT sparse representation
[x,v]=TQWT_SR_GMC_penalty_fun(y',Q_choose,r,J_fin,gamma,ka_choose,0);
y_GMC = itqwt(x,Q_choose,r,length(y));
figure,plot(y_GMC,'color',0.5*[1 1 1])
% xlim([0,10000])
% ylim([-3,3])
% figure,plot(t,sig1(1:length(y)),t,y_GMC,'r--')
% legend('故障特征信号','去噪信号')
% xlabel('时间[s]'),ylabel('幅值')

% figure,plot(t,sig1(1:length(y)),t,y_GMC,'r--')
% xlim([0.1,0.15])
% legend('故障特征成分','小波基函数')
% xlabel('时间[s]'),ylabel('幅值')
% set(gca,'xticklabel',[]),set(gca,'yticklabel',[])
% % figure,plot(t,y_GMC)
% % xlabel('时间[s]'),ylabel('幅值')
% % xlim([0,1])
% xlim([4000,6000])
% xlabel('采样时间')
% ylabel('幅值')
% set(gca,'xticklabel',[]),set(gca,'yticklabel',[])
% title('Q=1,r=3')
% 
% 
% figure,plot(sig_1)
% xlabel('时间[s]')
% ylabel('幅值')
% xlim([0,300])


blp=abs(fft(abs(hilbert(y_GMC))))/length(y_GMC)*2;
figure,plot(blp)
blp(1)=0;
% pl=(0:length(y_GMC)-1)/length(y_GMC)*fs;
% figure,plot(pl(1:round(length(y_GMC)/2)),blp(1:round(length(y_GMC)/2)))
% xlabel('频率 [Hz]'),ylabel('幅值')
xlim([0,350])


