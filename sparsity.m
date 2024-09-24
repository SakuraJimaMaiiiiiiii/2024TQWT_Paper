function sparsity = sparsity(x)
    sparsity=sqrt(sum(x.^2))/norm(x,1);
%      sparsity=sqrt(sum(x.^2))/sum(abs(x));
end

