
classdef Assembly_Truss < handle

    properties
        % Nodes
        node

        % Bars
        bar

    end

    methods
        % For a given input deformation find the global force vector and
        % the stiffness matrix
        [T,K]=SolveFK(obj,node,U)

        % Initialize the assembly
        % This will set currentU to be zero matrix
        % This will set theta_StressFree_Vec to be current theta value
        % This will set current external force to be zero vector
        InitializeAssembly(obj)

    end
end