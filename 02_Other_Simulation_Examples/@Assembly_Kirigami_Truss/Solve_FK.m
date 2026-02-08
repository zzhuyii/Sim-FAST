function [T,K]=Solve_FK(obj,U)

    [Tcst,Kcst]=obj.cst.Solve_FK(obj.node,U);
    [Trs,Krs]=obj.rot_spr_4N.Solve_FK(obj.node,U);
    [Tb,Kb]=obj.bar.Solve_FK(obj.node,U);

    T=Tcst+Trs+Tb;
    K=Kcst+Krs+Kb;

end