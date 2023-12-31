%% plot the deformation history of the simulated results.

function Plot_DeformedShape(obj,U)

View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;

assembly=obj.assembly;
undeformedNode=assembly.node.coordinates_Mat;

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

barConnect=assembly.bar.barConnect_Mat;
barNum=size(barConnect);
barNum=barNum(1);

for j=1:barNum
    node1=undeformedNode(barConnect(j,1),:);
    node2=undeformedNode(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color',[.7 .7 .7]);
end


wedgeConnect=assembly.wedge.wedgeConnect_Mat;
wedgeNum=size(wedgeConnect);
wedgeNum=wedgeNum(1);

for j=1:wedgeNum
    node1=undeformedNode(wedgeConnect(j,1),:);
    node2=undeformedNode(wedgeConnect(j,2),:);
    node3=undeformedNode(wedgeConnect(j,3),:);
    node4=undeformedNode(wedgeConnect(j,4),:);
    node5=undeformedNode(wedgeConnect(j,5),:);
    node6=undeformedNode(wedgeConnect(j,6),:);

    f=[1,2,3];
    v=[node1;node2;node3;];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

    f=[1,2,3];
    v=[node4;node5;node6;];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

    f=[1,2,3,4];
    v=[node1;node2;node5;node4];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

    f=[1,2,3,4];
    v=[node2;node3;node6;node5];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

    f=[1,2,3,4];
    v=[node1;node3;node6;node4];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)


end


deformNode=undeformedNode+U;

for j=1:barNum
    node1=deformNode(barConnect(j,1),:);
    node2=deformNode(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color','k');
end

for j=1:wedgeNum
    node1=deformNode(wedgeConnect(j,1),:);
    node2=deformNode(wedgeConnect(j,2),:);
    node3=deformNode(wedgeConnect(j,3),:);
    node4=deformNode(wedgeConnect(j,4),:);
    node5=deformNode(wedgeConnect(j,5),:);
    node6=deformNode(wedgeConnect(j,6),:);

    f=[1,2,3];
    v=[node1;node2;node3;];
    patch('Faces',f,'Vertices',v,'FaceColor','yellow')

    f=[1,2,3];
    v=[node4;node5;node6;];
    patch('Faces',f,'Vertices',v,'FaceColor','yellow')

    f=[1,2,3,4];
    v=[node1;node2;node5;node4];
    patch('Faces',f,'Vertices',v,'FaceColor','yellow')

    f=[1,2,3,4];
    v=[node2;node3;node6;node5];
    patch('Faces',f,'Vertices',v,'FaceColor','yellow')

    f=[1,2,3,4];
    v=[node1;node3;node6;node4];
    patch('Faces',f,'Vertices',v,'FaceColor','yellow')


end

