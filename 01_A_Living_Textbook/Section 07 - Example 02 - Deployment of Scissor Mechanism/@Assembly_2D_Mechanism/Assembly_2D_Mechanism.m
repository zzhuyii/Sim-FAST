
classdef Assembly_2D_Mechanism < handle

    properties
        % Nodes
        node

        % Bars
        bar

        % 3-node rot springs
        rot_spr_3N

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