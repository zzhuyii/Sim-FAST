function [passVec,modeStrList]=CheckAllBarDesignAISC(Nforce,A,E,Lc,r,Fy)
    nbar=size(A,1);
    passVec=zeros(nbar,1);
    modeStrList={};
    for i=1:nbar
    
        [pass,modeStr]=CheckBarDesignAISC(Nforce(i),A(i),E(i),Lc(i),r(i),Fy);
        passVec(i)=pass;
        modeStrList{i}=modeStr;
    end
end