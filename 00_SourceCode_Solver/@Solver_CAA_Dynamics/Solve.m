function Uhis=Solve(obj)

    % Input setup from loading controller
    assembly=obj.assembly;

    % support information
    supp=obj.supp;

    % time step
    dt=obj.dt;

    % external loading forces
    Fext=obj.Fext;

    % target folding angle
    rotSprTarget=obj.rotSprTargetAngle;

    % We need one more step for Fext 
    % The first step is with zeros
    A=size(Fext);
    step=A(1);
    nodeNum=A(2);

    % adjust the size of Fext 
    Fext0=zeros(1,nodeNum,3);
    Fext=cat(1,Fext0,Fext);

    % vector of every time step
    TimeVec=(1:step)*dt;

    
    % Set up storage         
    Uhis=zeros(step+1,nodeNum,3);
    VHis=Uhis;
    
    V0=zeros(size(assembly.node.current_U_mat));
    
    Uhis(1,:,:)=assembly.node.current_U_mat;
    VHis(1,:,:)=V0;
    
    % The static load that were previouslly applied
    currentAppliedForce=zeros(3*nodeNum,1);    
    for i=1:nodeNum
        currentAppliedForce(3*(i-1)+1:3*i) =...
            assembly.node.current_ext_force_mat(i,:);
    end  
   
    % Find the mass matrix of the system
    MassMat=assembly.node.FindMassMat();
    
    % Implement the explicit solver
    for i=1:step

        assembly.rot_spr_4N.theta_stress_free_vec=rotSprTarget(i,:)';
        [T,K]=assembly.Solve_FK(squeeze(Uhis(i,:,:)));

        [K,T]=Mod_K_For_Supp(K,supp,T);

        [K,Fexti]=Mod_K_For_Supp(K,supp,...
            reshape(squeeze(Fext(i,:,:))',[3*nodeNum,1]));
        [K,Fexti1]=Mod_K_For_Supp(K,supp,...
            reshape(squeeze(Fext(i+1,:,:))',[3*nodeNum,1]));
        
        [K,Vhisi]=Mod_K_For_Supp(K,supp,...
            reshape(squeeze(VHis(i,:,:))',[3*nodeNum,1]));
        [K,Uhisi]=Mod_K_For_Supp(K,supp,...
            reshape(squeeze(Uhis(i,:,:))',[3*nodeNum,1]));

        K=sparse(K);

            
        % Set up the damping matrix
        alpha=obj.alpha;
        beta=obj.beta;
        DampMat=alpha*MassMat+beta*K;
        
        
        % Solve the acceleration
        UDotDot_i=MassMat\(Fexti-DampMat*Vhisi-T);
        
        Kadjust=K+2/dt*DampMat+4/dt/dt*MassMat;
        dP_adjust=(Fexti1-Fexti)+2*DampMat*Vhisi...
            +MassMat*(4/dt*Vhisi+2*UDotDot_i);
        
        Uhisi1=Kadjust\dP_adjust+Uhisi;
        
        Vhisi1=2/dt*(Uhisi1-Uhisi)-Vhisi;
        
        Uhis(i+1,:,:)=reshape(Uhisi1,[3,nodeNum])';
        VHis(i+1,:,:)=reshape(Vhisi1,[3,nodeNum])';
        
        if rem(i,1000)==0
            fprintf('finish solving %d step \n',i);
        end       
        
    end

    Uhis=Uhis(1:step,:,:);

end