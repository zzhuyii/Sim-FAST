function [pass,modeStr,Pn,stressRatio]=Check_Truss_AISC(Ni,Ai,Ei,Lci,ri,Fy)

    if Ni>0
        % in tension 
        Pn_i=Fy*Ai;
        modeStr='Tension-Yield';

    else
        % in compression
        slender=Lci/ri;
        Fe=(pi^2*Ei)/(slender^2);

        % threshold for elastic/inelastic buckling
        lambda_lim=4.71*sqrt(Ei/Fy); 

        if slender<=lambda_lim
            Fcr=(0.658)^(Fy/Fe)*Fy;   % inelastic buckling
        else
            Fcr=0.877*Fe;             % elastic bukling
        end
        Pn_i=Fcr*Ai;
        modeStr='Compression-Buckling';
    end

    Pn=Pn_i;
    stressRatio=abs(Ni)/abs(Pn_i);
    pass=stressRatio<=1.0;

end