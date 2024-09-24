function [info_entropy, info_entropy_sum] = info_entropy(x)
%info_entropy为pij构成的参数矩阵
%info_entropy_sum为信息熵 即Ei
    info_entropy=zeros();
    info_entropy_sum=zeros();
 for i = 1:size(x,1)   
    for j =1:size(x,2)
        pij=x(i,j)/sum(x(i,:));
        info_entropy(i,j)=-pij*log(pij);
    end
%     info_entropy_sum(i)=sum(info_entropy(i,:));
 end
     for i=1:size(info_entropy,1)
        info_entropy_sum(i)=sum(info_entropy(i,:));
    end

end

  