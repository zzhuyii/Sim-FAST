%% Newton-Raphson Implicit Solver

classdef Solver_NR_Loading < handle
    properties

        % assembly of the structure
        assembly

        % storing the support information
        supp
        
        % the applied load
        load
        % a sampel code is the following
        % loadForce=3;
        % load=[29,0,0,loadForce;
        %       30,0,0,loadForce;
        %       31,0,0,loadForce;
        %       32,0,0,loadForce;];
        
        % the total number of incremental steps
        increStep=50
        
        % the tolerance for each iteration
        tol=1*10^-5
        
        % the maximum allowed iteration number
        iterMax=30        

        % the history of displacement field
        Uhis

        % factor for damped NR
        dampFactor=1

    end

    methods
        % Solve for the equilibrium results
        Uhis=Solve(obj);

    end
end