function Initialize(obj,node)

    rotSprIJKL=obj.node_ijkl_mat;
    numSpr=size(rotSprIJKL);
    numSpr=numSpr(1);

    % check the ijkl assignment to make sure they are the same direction
    for i=1:numSpr

        n1=rotSprIJKL(i,1);
        n2=rotSprIJKL(i,2);
        n3=rotSprIJKL(i,3);
        n4=rotSprIJKL(i,4);

        x1=node.coordinates_mat(n1,:);
        x2=node.coordinates_mat(n2,:);
        x3=node.coordinates_mat(n3,:);
        x4=node.coordinates_mat(n4,:);

        v1=x1-x2;
        v2=x3-x2;

        norm1=cross(v1,v2);

        if norm1(3)~=0
            if norm1(3)<0
                obj.node_ijkl_mat(i,2)=n3;
                obj.node_ijkl_mat(i,3)=n2;
            end

        elseif norm1(2)~=0
            if norm1(2)<0
                obj.node_ijkl_mat(i,2)=n3;
                obj.node_ijkl_mat(i,3)=n2;
            end
        else
            if norm1(1)<0
                obj.node_ijkl_mat(i,2)=n3;
                obj.node_ijkl_mat(i,3)=n2;
            end
        end

    end

    nodeSize=size(node.coordinates_mat());
    nodeSize=nodeSize(1);
    
    obj.theta_current_vec=obj.Solve_Theta(node,zeros(nodeSize,3));
    obj.theta_stress_free_vec=obj.theta_current_vec;

end

