%% 4-Node Rotational Spring Elements (directional stiffness)
% This rotational spring element is derived using analytical equations
% The rotational spring element is vectorized 
% The formulation is for a bi-linear elastic rotational spring
% This element uses different stiffness for different folding direction
% The rotational spring geometry is defined with 4 nodes

classdef Vec_Elements_RotSprings_4N_Directional_Smooth < handle

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

        % the MV assignment of each element (Ns*1)
        % 0 means mountain fold: fold downward
        % 1 means valley fold: fold upward
        mv_vec 

        % MV stiffness (Ns*1)
        % Here we assume that the folding stiffness is directional, the
        % element is more likely to fold in one direction compare to the
        % other direction. This is used to consider contact related
        % behaviors for thick origami hinges
        mv_factor_vec

        % Threshold for local penertration prevention
        % The theta1 and theta2 are terms for initiating the penertration 
        % prevention. When folding angle is smaller than theta1 or bigger 
        % than theta 2, an additional forces is added. Please check the 
        % Liu and Paulino RSPA paper for details.
        theta1=0.1*pi; 
        theta2=2*pi-0.1*pi;

    end

    methods
        % This function initialize the rotational spring element
        Initialize(obj,node)

        % Calculate the rotational angle of rotational springs
        [theta]=Solve_Theta(obj,node,U);

        % Calculate the local moment and stiffness of the rotational spring
        % elements. This function gives the constitutive relationship of 
        % rotational spring elements
        [Mspr,Cspr]=Solve_Moment(obj,theta)

        % Calculate the global internal force of the spring element
        [Tspr]=Solve_Global_Force(obj,node,U,Mspr)

        % Calculate the gloabl stiffness matrix of the spring element
        [Kspr]=Solve_Global_Stiff(obj,node,U,Mspr,Cspr)        

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the bar elements (making use of the above four 
        % functions). 
        [Tspr,Kspr]=Solve_FK(obj,node,U)


    end
end
