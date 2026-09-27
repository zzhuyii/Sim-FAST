function [Tspr,Kspr]=Solve_FK(obj,node,U)
    
        [theta]=Solve_Theta(obj,node,U);
        [Mspr,Cspr]=Solve_Moment(obj,theta);
        [Tspr]=Solve_Global_Force(obj,node,U,Mspr);
        [Kspr]=Solve_Global_Stiff(obj,node,U,Mspr,Cspr);

        obj.theta_current_vec=theta;
        obj.energy_current_vec=1/2*obj.rot_spr_K_vec.*...
            (theta-obj.theta_stress_free_vec).*...
            (theta-obj.theta_stress_free_vec);

end