function [Fj,J] = chooseband(x)
%该函数输出融合的特征参数Fj及挑选后的子带J
    xij_norm=zeros();
    Wi= weight(x);
%归一化参数矩阵    
for i = 1:size(x,1)   
    for j =1:size(x,2)
        xij_norm(i,j)=(x(i,j)-min(x(i,:)))/(max(x(i,:))-min(x(i,:)));
    end
end
    Fj=zeros();
%获得融合特征参数
    for j=1:size(x,2)
        Fj(j)=Wi*xij_norm(:,j);
    end
%挑选最大值对应的子带 即J的大小
    [~,J]=max(Fj);  %~为矩阵最大值，J为最大值的下表 由于main中范围为1:1:10 则返回的位置即为Q值
end

