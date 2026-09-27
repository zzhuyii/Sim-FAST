function  [Mspr,Cspr]= Solve_Moment(obj,theta)

    theta1=obj.theta1;
    theta2=obj.theta2;

    sprRotK=obj.rot_spr_K_vec;
    spr_StressFree=obj.theta_stress_free_vec;

    Mspr=zeros(size(theta));
    Cspr=zeros(size(theta));

    A=size(theta);
    N=A(1);
   
    % These codes are not vectorized yet
    % The constitutive relationship is not vectorzied because the 
    % if loop is not easily vectorized. Moreover, the computation is 
    % not too time consuming.
    for i=1:N           
        
        % if theta(i)<pi && obj.mv_vec(i)==0
        %     tempK = sprRotK(i);
        % elseif theta(i)>pi && obj.mv_vec(i)==1
        %     tempK = sprRotK(i);
        % else
        %     tempK = sprRotK(i)*obj.mv_factor_vec(i);
        % end
        
        if obj.mv_vec(i)==0
            % theta > Pi 
            factor=(tanh((pi-theta(i))*20)+1)*obj.mv_factor_vec(i)+1;
            tempK = sprRotK(i)*factor;
        else 
            % theta < Pi 
            factor=(tanh((theta(i)-pi)*20)+1)*obj.mv_factor_vec(i)+1;
            tempK = sprRotK(i)*factor;
        end

        if theta(i)<theta1
            ratio=(sec(pi*(theta(i)-theta1)/2/theta1))^2;
            Cspr(i)=tempK*ratio;
            Mspr(i)=tempK*(theta(i)-spr_StressFree(i))...
                +(2*theta1*tempK/pi)*...
                tan(pi*(theta(i)-theta1)/2/theta1)-...
                tempK*theta(i)+tempK*theta1;
     
        elseif theta(i)>theta2
            ratio=(sec(pi*(theta(i)-theta2)/(4*pi-2*theta2)))^2;
            Cspr(i)=tempK*ratio;
            Mspr(i)=tempK*(theta(i)-spr_StressFree(i))...
                +(2*(2*pi-theta2)*tempK/pi)...
                *tan(pi*(theta(i)-theta2)/(4*pi-2*theta2))-...
                tempK*theta(i)+tempK*theta2;
        else   
            Cspr(i)=tempK; 
            Mspr(i)=tempK*(theta(i)-spr_StressFree(i));
        end

    end
        
end
