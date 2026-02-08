
classdef Assembly_Scissor_Bridge < handle

    properties
        % Nodes
        node

        % CST Elements
        cst

        % Bar Elements
        bar

        % 3 Node Rotational Springs
        rot_spr_3N

        % 4 Node Rotational Springs
        rot_spr_4N

        % Act Bar
        actBar
        
    end
    
    methods

        % For a given input deformation find the global force vector and
        % the stiffness matrix
        [T,K]=Solve_FK(obj,node,U)

        % Initialize the assembly
        % This will set currentU to be zero matrix
        % This will set theta_StressFree_Vec to be current theta value
        % This will set current external force to be zero vector
        Initialize_Assembly(obj)

    end
end