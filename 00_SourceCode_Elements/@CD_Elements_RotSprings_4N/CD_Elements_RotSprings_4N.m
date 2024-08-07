%% 4-Node Rotational Spring Elements
% This rotational spring element is derived using central difference method
% The formulation is for a linear-elastic rotational spring
% The rotational spring geometry is defined with 4 nodes

classdef CD_Elements_RotSprings_4N < handle

    properties
        % Node number of the four node used to define the rotational 
        % spring elements (Ns*4)
        node_ijkl_mat

        % Rotational stiffness of each element (Ns*1)
        rot_spr_K_vec  

        % Stress-free angle of the spring
        theta_stress_free_vec;

        % Current Theta
        theta_current_vec;

        % Current Strain Energy
        energy_current_vec;

        % step for central difference
        delta=10^-8;
        
    end

    methods
        % This function initialize the springs
        Initialize(obj,node)

        % Calculate the potential of the rotational spring
        PE=Potential(obj,theta,theta0,K);

        % Calculate the rotational angle of rotational springs
        theta=Solve_Theta(obj,X);

        % Calculate the local force vector of the rotational spring
        [Flocal]=Solve_Local_Force(obj,X,theta0,K)

        % Calculate the global internal force of all spring elements
        [Trs]=Solve_Global_Force(obj,node,U)

        % Calculate the local stiffness of the rotational spring
        [Klocal]=Solve_Local_Stiff(obj,X,theta0,K)

        % Calculate the gloabl stiffness matrix of the spring element
        [Krs]=Solve_Global_Stiff(obj,node,U)        

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the rotational spring elements (making use of the 
        % above four functions). 
        [Trs,Krs]=Solve_FK(obj,node,U)

    end
end
