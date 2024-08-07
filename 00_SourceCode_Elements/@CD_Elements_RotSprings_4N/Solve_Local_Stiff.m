function [Klocal]=Solve_Local_Stiff(obj,X,theta0,K)

    % The rotational spring elemement has 4 nodes, the hessian should be 
    % a 12 by 12 matrix. Local stiffness is the numerical hessian of the 
    % potential.

    Klocal=zeros(12,12);
    
    % We use this small step delta to calculate numerical Hessian
    delta=obj.delta;

    % for each element in the matrix, we do the following calculation    
    for i=1:12

        % find the index
        if mod(i,3)==0
            index1=3;
            index2=i/3;
        else
            index1=mod(i,3);
            index2=floor(i/3)+1;
        end

        % forward step, where the coordinate is moved in the selected
        % direction by +delta
        tempXfor=X;
        tempXfor(index2,index1)=tempXfor(index2,index1)+delta;
        
        % backward step, where the coordinates is moved in the selected
        % direction by -delta
        tempXback=X;
        tempXback(index2,index1)=tempXback(index2,index1)-delta;
    
        % This is the central deference equation. Please note that we can 
        % make use of the force function for this calculation. 
        Klocal(i,:)=1/2/delta*(obj.Solve_Local_Force(tempXfor,theta0,K)-...
            obj.Solve_Local_Force(tempXback,theta0,K));


    end    
end