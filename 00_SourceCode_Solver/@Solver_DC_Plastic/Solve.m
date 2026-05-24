%% Displacement controled method method. 

function [Uhis,Fhis]=Solve(obj)

    % initialize and set up storage matrix for Uhis
    increStep=obj.increStep;
    tol=obj.tol;
    iterMax=obj.iterMax;
    lambdaBar=obj.lambdaBar;  
    selectedRefDisp=obj.selectedRefDisp;
    
    supp=obj.supp;
    load=obj.load;   
  
    assembly=obj.assembly;
    U=assembly.node.current_U_mat;

    A=size(U);
    NodeNum=A(1);
    Uhis=zeros(increStep,NodeNum,3);
    Fhis=zeros(increStep,NodeNum*3);

    % Find the external forces that is currently applied on the structure
    currentAppliedForce=zeros(3*NodeNum,1);    
    for i=1:NodeNum
        currentAppliedForce(3*(i-1)+1:3*i) =...
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
     
    fprintf('Loading Analysis Start');
 
    pload=loadVec;
    up1=zeros(NodeNum*3,increStep);
    lambda=1;

    selectedReferenceNodeNum=selectedRefDisp(1)*3-3+selectedRefDisp(2);
    
    % This code can automatically adjust the loading step
    for i=1:increStep

        step=1;     
        fprintf('Icrement = %d\n',i);

        % find the internal force and stiffness of system
        [T,K]=assembly.Solve_FK(U);

        % calculate the unbalanced force
        unbalance=currentAppliedForce+lambda*loadVec-T; 

        [K,unbalance]=Mod_K_For_Supp(K,supp,unbalance);
        K=sparse(K);                         

        up1(:,i)=K\loadVec;  
        
        % If there is no deformation, we need GSP to be zero
        if norm(up1)==0
            GSP=1;
            sig=1;
        else
            if i==1
                GSP=1;
                sig=1;
            else
                GSP=(up1(:,1)'*up1(:,1))/(up1(:,i)'*up1(:,i));
                sig=sign(up1(:,i-1)'*up1(:,i))*sig;        
            end
        end
        
        % we can deactive the GSP scalling
        % GSP=1;
        % sig=1;
        
        dLambda=sig*lambdaBar*sqrt(abs(GSP));
        lambda=lambda+dLambda;
        pload=pload+dLambda*loadVec;    
        dUtemp=dLambda*up1(:,i);

        for j=1:NodeNum
            U((j),:)=U((j),:)+dUtemp(3*j-2:3*j)';
        end  
        R=norm(dUtemp);
        fprintf('    Iteration = %d, R = %e\n',step,R);
        
        step=step+1;       

        while and(step<iterMax,R>tol)

            % find the internal force and stiffness of system
            [T,K]=assembly.Solve_FK(U);
    
            % calculate the unbalanced force
            unbalance=currentAppliedForce+lambda*loadVec-T; 
    
            [K,unbalance]=Mod_K_For_Supp(K,supp,unbalance);
            K=sparse(K);  

            up=K\loadVec;
            ur=K\unbalance;

            dLambda=-ur(selectedReferenceNodeNum)/up(selectedReferenceNodeNum);

            lambda=lambda+dLambda;
            pload=pload+dLambda*loadVec;
            
            dUtemp=dLambda*up+ur;
            for j=1:NodeNum
                U((j),:)=U((j),:)+dUtemp(3*j-2:3*j)';
            end

            R=norm(dUtemp);
            fprintf('    Iteration = %d, R = %e\n',step,R);
         
            step=step+1;  
        end 

        % Store the found equilibrium
        Uhis(i,:,:)=U;
        Fhis(i,:)=T;
   
    end    
end

