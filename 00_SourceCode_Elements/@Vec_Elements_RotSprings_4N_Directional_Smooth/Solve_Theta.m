function [theta]=Solve_Theta(obj,node,U)

    nodalCoordinate=node.coordinates_mat;
    sprIJKL=obj.node_ijkl_mat;

    spr_i=sprIJKL(:,1);
    spr_j=sprIJKL(:,2);
    spr_k=sprIJKL(:,3);
    spr_l=sprIJKL(:,4);
    
    nodei=nodalCoordinate(spr_i,:)+U(spr_i,:);
    nodej=nodalCoordinate(spr_j,:)+U(spr_j,:);
    nodek=nodalCoordinate(spr_k,:)+U(spr_k,:);
    nodel=nodalCoordinate(spr_l,:)+U(spr_l,:);
    
    rij=nodei-nodej;
    rkj=nodek-nodej;
    rkl=nodek-nodel;
    
    m=cross(rij,rkj,2);
    n=cross(rkj,rkl,2);
    
    d_m_kl=dot(m,rkl,2);
    zero_index=find(d_m_kl==0);

    yita=sign(d_m_kl);
    yita(zero_index)=1;
    
    m_norm=sqrt(dot(m,m,2));
    n_norm=sqrt(dot(n,n,2));
    
    theta=mod(yita.*real(acos(dot(m,n,2)./m_norm./n_norm)),2*pi);
    
end