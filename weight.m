function Wi= weight(x)
     %x为参数矩阵 y为矩阵长度
     %输出权重Wi Wi为长度为参数数目的行向量
     [~, Ei] = info_entropy(x);
     Wi=zeros();
     for i=1:size(x,1)
        Wi(i)=1-Ei(i)/log(size(x,2));
     end
end

