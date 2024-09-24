function tqwt_eng=calculate_energy(Q,r,J)
    J1=1;
    J2=J;
    N=2000;
    beta = 2/(Q+1);
    alpha = 1-beta/r;
    wlets = ComputeWavelets(N,Q,r,J2,'radix2');
    for i=1:J+1
        tqwt_eng(i)=sum(abs(wlets{i}).^2);
    end


end
