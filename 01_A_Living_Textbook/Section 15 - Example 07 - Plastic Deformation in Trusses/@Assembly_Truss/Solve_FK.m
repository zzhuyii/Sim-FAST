%% Find the global internal force and stiffness matrix of the assembly

function [T,K]=Solve_FK(obj,U)

    [Tbar,Kbar]=obj.bar.Solve_FK(obj.node,U);
    %[Tabar,Kabar]=obj.actBar.Solve_FK(obj.node,U);
    
    % T=Tbar+Tabar;
    % K=Kbar+Kabar;

    T=Tbar;
    K=Kbar;
end