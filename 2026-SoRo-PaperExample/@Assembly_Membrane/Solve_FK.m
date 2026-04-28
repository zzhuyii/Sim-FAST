function [T,K]=Solve_FK(obj,U)

    [Tcst,Kcst]=obj.cst.Solve_FK(obj.node,U);
    [Trs,Krs]=obj.rotSpr.Solve_FK(obj.node,U);
    % [Tab,Kab]=obj.actBar.Solve_FK(obj.node,U);

    %T=Tcst+Trs+Tab;
    %K=Kcst+Krs+Kab;

    T=Tcst+Trs;
    K=Kcst+Krs;
    
end