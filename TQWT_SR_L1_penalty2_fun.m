function x1=TQWT_SR_L1_penalty2_fun(y,Q,r,J,gamma,ka,state)
    tqwt_eng=calculate_energy(Q,r,J);
    max_iter=4000; % maximum of the iteration
    tol=1e-3;
    soft = @(x, T) max(1 - T./abs(x), 0) .* x;
    mu = 1.9 / ( 1 * max( 1,  gamma / (1-gamma) ) );
    
    x=tqwt(zeros(size(y)),Q,r,J);
   
   
    iter=0;
    del_x=inf;
    old_x=x{1}/tqwt_eng(1);
    for i=1:J
        old_x=[old_x(:);x{i+1}(:)/tqwt_eng(i+1)];
    end
    x1=cell(J+1,1);
     
        for i=1:J+1
          x1{i}=x{i};
        end
        
    while (iter<max_iter)&&(del_x>tol)
        iter=iter+1
        
        Ax1=itqwt(x1,Q,r,length(y));
        ATAx1=tqwt(Ax1,Q,r,J);
        ATy=tqwt(y,Q,r,J);
        for i=1:J+1
            w{i}=x1{i}-mu*ATAx1{i}+mu*ATy{i};
            
        end
        
        clear norm_w
        norm_w=w{1}/tqwt_eng(1);
        for i=1:J
            norm_w=[norm_w(:);w{i+1}(:)/tqwt_eng(i+1)];
        end
        if state==0   % method used in the literature
            k=50;
            n_w=sort(norm_w,1,'descend');
            T=n_w(k);
        else
            norm_1_w=abs(norm_w);
            T=ka*(max(norm_1_w)-min(norm_1_w))+min(norm_1_w);
        end
        
        for i=1:J+1
            x1{i}=soft(w{i},T*tqwt_eng(i));
           
        end
        del_x=max(abs(norm_w-old_x))/max(abs(old_x));
        old_x=norm_w;
        
    end
end