function [pass,modeStr]=CheckBarDesignAISC(Ni,Ai,Ei,Lci,ri,Fy)
    
    tiny=1e-12;
    
    if Ni>0
        % in tension 
        Pn_i=Fy*Ai;
        modeStr='Tension-Yield';
        slender=NaN;  
        Fe=NaN; 
        Fcr=NaN;   % don't calculate buckling in tension
    else
        % in compression
        slender=Lci/ri;
        Fe=(pi^2*Ei)/(slender^2);
        lambda_lim=4.71*sqrt(Ei/Fy); % threshold for elastic/inelastic buckling

        if slender<=lambda_lim
            Fcr=(0.658)^(Fy/Fe)*Fy;   % inelastic buckling
        else
            Fcr=0.877*Fe;             % elastic bukling
        end
        Pn_i=Fcr*Ai;
        modeStr='Compression-Buckling';
    end

    Pn=Pn_i;
    util=abs(Ni)/max(Pn_i,tiny);
    pass=util<=1.0;

end