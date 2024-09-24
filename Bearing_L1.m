%% import signal
clear
clc
close all
%% chosen signal
% original signal
M = csvread('60.csv',1,0);
sig0=M(:,1);
figure,plot(sig0,'color',0.5*[1 1 1])
xlim([0,10000])
ylim([-9,9])
fs=26500;  % sampling freq.

% 20000length signal
y = sig0;
y = y(1:20000);
figure,plot(y)
xlabel('时间 [s]'),ylabel('幅值')
xlim([1*10^4,1.5*10^4])
xticklabels({'1','1.05','1.1','1.15','1.2','1.25','1.3','1.35','1.4','1.45','1.5'})

%% TQWT parameters
% High Q-factor wavelet transform
gamma=0.85;
%%% w = tqwt(x,Q,r,J);              % TQWT
%%%y = itqwt(w,Q,r,length(x));             % Inverse TQWT




%% Q,R inter
Q_range = 1:1:10;               
r_range = 3:1:5;
ka_range = 0.3:0.1:1;
N=200;
ESHE = zeros();
T1 = 107.9;
f=1/T1;
delt_p=1.5;
fz1 = zeros();
ff1 = zeros();
for i = 1 :length(Q_range)
    Q= Q_range(i);
    for j = 1:length(r_range)
        r = r_range(j);
        J = min(10,computeJmax(N,Q,r));
        for k=1:length(ka_range)
            ka = ka_range(k);
            x=TQWT_SR_L1_penalty2_fun(y,Q,r,J,gamma,ka,0);
            y_GMC = itqwt(x,Q,r,length(y));
%%%---%%% envelope spectrum harmonic energy
            blp=abs(fft(abs(hilbert(y_GMC))))/length(y_GMC)*2;
            ff0=(0:length(y_GMC)-1)/length(y_GMC)*fs;
            ff=ff0(1:round(length(y_GMC)/2));
            bl=blp(1:round(length(y_GMC)/2));
            bl(1)=0;
            jg=fix(f/fs*length(1:round(length(y_GMC)/2)));
            [pks,loc]=findpeaks(bl,'minpeakdistance',round(jg/10));
            f_temp=f;
            clear fz1 ff1
            for ij=1:5
                [~,zb0]=find(ff>=f_temp-delt_p,1,'first');
                [~,yb0]=find(ff<=f_temp+delt_p,1,'last');
                [~,loc1]=find(loc>=zb0 & loc<=yb0);
                [~,loc10]=min(abs(ff(loc(loc1))-f_temp));
                if ~isempty(loc10)
                    fz1(ij)=bl(loc(loc1(loc10)));
                    ff1(ij)=ff(loc(loc1(loc10)));
                    f_temp=f+ff(loc(loc1(loc10)));
                else
                    fz1(ij)=0;
                    f_temp=f+f_temp;
                end
            end
            ESHE(i,j,k)=sum(fz1(1:3).^2);      
        end
    end
end

%% parm sort
c = reshape(ESHE,1,[]);
c =sort(c,'descend');
figure,plot(c) 
ylabel('ESHE')
xlabel('数据长度')


%% choose Q,r,J
[i,j]=max(ESHE);
i=squeeze(i);
j=squeeze(j);
[i_,j_]=max(i);
[i__,j__]=max(i_);
max_k=j__;
max_j=j_(max_k);
max_i=j(max_j,max_k);
Q = Q_range(max_i);
r = r_range(max_j);
ka = ka_range(max_k);
J = min(10,computeJmax(N,Q,r));


%% TQWT sparse representation
x=TQWT_SR_L1_penalty2_fun(y,Q,r,J,gamma,ka,0)
y_GMC = itqwt(x,Q,r,length(y));
figure,plot(y_GMC,'color',0.5*[1 1 1]);
legend('故障特征信号','去噪信号')
xlabel('时间[s]'),ylabel('幅值')
xlim([0,0.4])
figure,plot(y_GMC,'color',0.5*[1 1 1])
xlabel('时间[s]'),ylabel('幅值')
blp=abs(fft(abs(hilbert(y_GMC))))/length(y_GMC)*2;
figure,plot(blp,'color',0.5*[1 1 1])
xlim([0,450]);
ylim([-2.2]);
blp(1)=0;
pl=(0:length(y_GMC)-1)/length(y_GMC)*fs;
figure,plot(pl(1:round(length(y_GMC)/2)),blp(1:round(length(y_GMC)/2)),'color',0.5*[1 1 1])
xlabel('频率 [Hz]'),ylabel('幅值')
xlim([0,450])

set(gca,'ytick',[],'yticklabel',[])
