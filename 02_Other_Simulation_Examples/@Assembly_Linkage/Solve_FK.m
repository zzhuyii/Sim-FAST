function [T,K]=Solve_FK(obj,U)

    [Tbar,Kbar]=obj.bar.Solve_FK(obj.node,U);
    [Trs,Krs]=obj.rot_spr_4N.Solve_FK(obj.node,U);
    
    T=Tbar+Trs;
    K=Kbar+Krs;
end