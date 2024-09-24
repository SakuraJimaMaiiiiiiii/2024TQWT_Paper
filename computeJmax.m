function  Jmax = computeJmax(N,Q,r)
          Jmaxdouble = [log(N/4/(Q+1))]/[log((Q+1)/(Q+1-2/r))];
          Jmax = floor(Jmaxdouble);
end
