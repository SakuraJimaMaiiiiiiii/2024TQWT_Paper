function ESHE = ES_EHNR(y)
%%%---%%% envelope spectrum harmonic energy
            fz1 = zeros();
            ff1 = zeros();
            fs=26500;
            f=1/107.9;
            delt_p=1.5;
            blp=abs(fft(abs(hilbert(y))))/length(y)*2;
            ff0=(0:length(y)-1)/length(y)*fs;
            ff=ff0(1:round(length(y)/2));
            bl=blp(1:round(length(y)/2));
            bl(1)=0;
            jg=fix(f/fs*length(1:round(length(y)/2)));
            [~,loc]=findpeaks(bl,'minpeakdistance',round(jg/10));
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
            ESHE=sum(fz1(1:3).^2);
end

