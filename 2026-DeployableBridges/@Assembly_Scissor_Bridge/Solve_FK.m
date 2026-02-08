function [T, K] = Solve_FK(obj, U)  

    [Tcst, Kcst] = obj.cst.Solve_FK(obj.node, U);
    [Tbar, Kbar] = obj.bar.Solve_FK(obj.node, U);
    [Tspr3, Kspr3] = obj.rot_spr_3N.Solve_FK(obj.node, U);
    [Tspr4, Kspr4] = obj.rot_spr_4N.Solve_FK(obj.node, U);
    [Tab, Kab] = obj.actBar.Solve_FK(obj.node, U);

    T = Tcst + Tbar + Tspr3 + Tspr4 + Tab;
    K = Kcst + Kbar + Kspr3 + Kspr4 + Kab;
end