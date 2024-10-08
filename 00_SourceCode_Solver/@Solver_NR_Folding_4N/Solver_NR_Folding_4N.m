%% Newton-Raphson Implicit Solver

classdef Solver_NR_Folding_4N  < handle
    properties

        % assembly of the structure
        assembly

        % storing the support information
        supp
               
        % the total number of incremental steps
        increStep=50
        
        % the tolerance for each iteration
        tol=1*10^-5
        
        % the maximum allowed iteration number
        iterMax=30
     
        % target stress free angle of creases
        targetRot
        
        % the history of displacement field
        Uhis
        
    end    

    methods
        % Solve for the equilibrium position
        Uhis=Solve(obj);

    end
end