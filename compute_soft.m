
function out = compute_soft(x,thres,alpha)
        out=zeros(size(x));
        ind1=find(abs(x)<=alpha*thres);
        ind2=find(abs(x)>alpha*thres & abs(x)<thres/alpha);
        ind3=find(abs(x)>=thres/alpha);
        
        
        out(ind1)=0;
        out(ind2)=x(ind2).*(abs(x(ind2))-alpha*thres)/(thres/alpha-thres*alpha);
        out(ind3)=x(ind3);


         
