%% Find the global internal force and stiffness matrix of the assembly

function [T,K]=Solve_FK(obj,U)

    [Tbar,Kbar]=obj.bar.Solve_FK(obj.node,U);
    [Tabar,Kabar]=obj.actBar.Solve_FK(obj.node,U);
    [Tspr,Kspr]=obj.rot_spr_4N.Solve_FK(obj.node,U);

    [Tcst,Kcst]=obj.cst.Solve_FK(obj.node,U);
    [Trs,Krs]=obj.rot_spr_4N.Solve_FK(obj.node,U);
    
    T=Tbar+Tabar+Tspr+Tcst+Trs;
    K=Kbar+Kabar+Kspr+Kcst+Krs;
end