%% plot the deformation history of the simulated results.

function Plot_DeformedShape(obj,Uundeformed,Udeformed)

View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;

assembly=obj.assembly;
undeformedNode=assembly.node.coordinates_Mat+Uundeformed;

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

deformNode=assembly.node.coordinates_Mat+Udeformed;

for j=1:barNum
    node1=deformNode(barConnect(j,1),:);
    node2=deformNode(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color','k','Linewidth',1.4);
end

panelNum=size(obj.panelConnection);
panelNum=panelNum(2);

for i=1:panelNum
    nodeNumVec=obj.panelConnection{i};

    f=[];
    v=[];
    for j=1:length(nodeNumVec)
        f=[f,j];
        v=[v;deformNode(nodeNumVec(j),:)];
    end

    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',0.2)

end

