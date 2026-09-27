function [Tspr]=Solve_Global_Force(obj,node,U,M)
    
    A=size(U);
    NodeNum=A(1);    
    Tspr=zeros(3*NodeNum,1);

    nodalCoordinates=node.coordinates_mat;
    sprIJKL=obj.node_ijkl_mat;

    spr_i=sprIJKL(:,1);
    spr_j=sprIJKL(:,2);
    spr_k=sprIJKL(:,3);
    spr_l=sprIJKL(:,4);

    spr_i=nonzeros(spr_i);
    spr_j=nonzeros(spr_j);
    spr_k=nonzeros(spr_k);
    spr_l=nonzeros(spr_l);

    nodei=nodalCoordinates(spr_i,:)+U(spr_i,:);
    nodej=nodalCoordinates(spr_j,:)+U(spr_j,:);
    nodek=nodalCoordinates(spr_k,:)+U(spr_k,:);
    nodel=nodalCoordinates(spr_l,:)+U(spr_l,:);

    rij=(nodei-nodej);
    rkj=(nodek-nodej);
    rkl=(nodek-nodel);

    m=cross(rij,rkj,2);
    n=cross(rkj,rkl,2);  

    m_square=dot(m',m')';
    n_square=dot(n',n')';

    rkj_square=dot(rkj',rkj')';
    rkj_norm=sqrt(rkj_square);

    parti=(rkj_norm./m_square).*m;
    partl=-rkj_norm./n_square.*n;
    partj=((dot(rij',rkj')'./rkj_norm./rkj_norm-1)).*...
        parti-dot(rkl',rkj')'./rkj_norm./rkj_norm.*partl;
    partk=((dot(rkl',rkj')'./rkj_norm./rkj_norm-1)).*...
        partl-dot(rij',rkj')'./rkj_norm./rkj_norm.*parti;

    gradient=[parti,partj,partk,partl];
    
    M_temp=M;    
    localT=M_temp.*gradient;        
    
    index1=3*spr_i-2;
    index2=3*spr_j-2;
    index3=3*spr_k-2;
    index4=3*spr_l-2;

    index=[index1, index1+1, index1+2,...
             index2, index2+1, index2+2,...
             index3, index3+1, index3+2,...
             index4, index4+1, index4+2];
         
    
    index=index(:);
    localT=localT(:);
    
    for i=1:length(index)
        Tspr(index(i))=Tspr(index(i))+localT(i);        
    end
end
