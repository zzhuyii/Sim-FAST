function [T,K]=Solve_FK(obj,U)

    [Tcst,Kcst]=obj.cst.Solve_FK(obj.node,U);
    [Trs_D,Krs_D]=obj.rot_spr_4N_D.Solve_FK(obj.node,U);

    if isempty(obj.rot_spr_4N)
    else
        [Trs,Krs]=obj.rot_spr_4N.Solve_FK(obj.node,U);
        [Tzlspr,Kzlspr]=obj.zlspr.Solve_FK(obj.node,U);
    end

    if isempty(obj.rot_spr_4N)
        T=Tcst+Trs_D;
        K=Kcst+Krs_D;
    else
        T=Tcst+Trs+Tzlspr+Trs_D;
        K=Kcst+Krs+Kzlspr+Krs_D;
    end

end