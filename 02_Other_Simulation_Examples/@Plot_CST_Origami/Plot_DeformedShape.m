%% plot the deformation history of the simulated results.

function Plot_DeformedShape(obj,U)

View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;

assembly=obj.assembly;
undeformedNode=assembly.node.coordinates_mat;

figure;

set(gcf, 'color', 'white');
set(gcf,'position',[obj.x0,obj.y0,obj.width,obj.height])

view(View1,View2); 
set(gca,'DataAspectRatio',[1 1 1])

A=size(Vsize);
if A(1)==1    
    axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
else
    axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
end

node0=assembly.node.coordinates_mat;
A=size(assembly.node.coordinates_mat);

cstIJK=obj.assembly.cst.node_ijk_mat;
panelNum=size(cstIJK);
panelNum=panelNum(1);

for k=1:panelNum
    nodeNumVec=cstIJK(k,:);

    f=[];
    v=[];
    for j=1:length(nodeNumVec)
        f=[f,j];
        v=[v;node0(nodeNumVec(j),:)];
    end

    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',0.2)

end

deformNode=undeformedNode+U;

for k=1:panelNum
    nodeNumVec=cstIJK(k,:);

    f=[];
    v=[];
    for j=1:length(nodeNumVec)
        f=[f,j];
        v=[v;deformNode(nodeNumVec(j),:)];
    end

    patch('Faces',f,'Vertices',v,'FaceColor','yellow','FaceAlpha',0.2)

end


