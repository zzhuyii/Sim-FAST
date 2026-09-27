function [T,K]=Solve_FK(obj,U)

    [Tcst,Kcst]=obj.cst.Solve_FK(obj.node,U);
    [Trs,Krs]=obj.rotSpr.Solve_FK(obj.node,U);
    
    T=Tcst+Trs;
    K=Kcst+Krs;
end