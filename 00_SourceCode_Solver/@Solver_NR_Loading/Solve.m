%% Newton-Raphson solver for loading

function [Uhis]=Solve(obj)

    % initialize and set up storage matrix for Uhis
    increStep=obj.increStep;
    tol=obj.tol;
    iterMax=obj.iterMax;
    
    supp=obj.supp;
    load=obj.load;   
  
    assembly=obj.assembly;
    U=assembly.node.current_U_mat;

    A=size(U);
    NodeNum=A(1);
    Uhis=zeros(increStep,NodeNum,3);

    fprintf('Loading Analysis Start');
    
    % Find the external forces that is currently applied on the structure
    currentAppliedForce=zeros(3*NodeNum,1);    
    for i=1:NodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = ...
            assembly.node.current_ext_force_mat(i,:);
    end  
    
    % Assemble the load vector
    A=size(load);
    loadSize=A(1);
    loadVec=zeros(3*NodeNum,1);
    for i=1:loadSize
        TempNodeNum=(load(i,1));
        loadVec(TempNodeNum*3-2)=load(i,2);
        loadVec(TempNodeNum*3-1)=load(i,3);
        loadVec(TempNodeNum*3-0)=load(i,4);
    end

    for i=1:increStep
        step=1;
        R=1;
        lambda=i;
        fprintf('Icrement = %d\n',i);
        
        
        while and(step<iterMax,R>tol)

            % find the internal force and stiffness of system
            [T,K]=assembly.Solve_FK(U);

            % calculate the unbalanced force
            unbalance=currentAppliedForce+lambda*loadVec-T; 
            
            [K,unbalance]=Mod_K_For_Supp(K,supp,unbalance);
            K=sparse(K);
                        
            dUtemp=(K\unbalance);
            for j=1:NodeNum
                U((j),:)=U((j),:)+dUtemp(3*j-2:3*j)'*obj.dampFactor;
            end
            R=norm(dUtemp);
            fprintf('    Iteration = %d, R = %e\n',step,R);
            step=step+1;        
        end

        Uhis(i,:,:)=U;
    end  

    assembly.node.current_U_mat=U;
end