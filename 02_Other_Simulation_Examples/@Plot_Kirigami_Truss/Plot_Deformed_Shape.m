%% plot the deformation history of the simulated results.

function Plot_Deformed_Shape(obj,U)

% Obtain Information
View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;
assembly=obj.assembly;
undeformedNode=assembly.node.coordinates_mat;

%% Set up the graphic space
figure
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

%% Plot the undeformed cst element 
node0=assembly.node.coordinates_mat;
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

%% Plot the undeformed bar element
barConnect=obj.assembly.bar.node_ij_mat;
barNum=size(barConnect);
barNum=barNum(1);

for j=1:barNum
    node1=node0(barConnect(j,1),:);
    node2=node0(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color',[0.5,0.5,0.5] );
end

%% Plot the deformed CST element
deformNode=undeformedNode+U;

for k=1:panelNum
    nodeNumVec=cstIJK(k,:);
    f=[];
    v=[];
    for j=1:length(nodeNumVec)
        f=[f,j];
        v=[v;deformNode(nodeNumVec(j),:)];
    end
    patch('Faces',f,'Vertices',v,'FaceColor','yellow')

end

%% Plot the deformed bar element
barConnect=obj.assembly.bar.node_ij_mat;
barNum=size(barConnect);
barNum=barNum(1);

for j=1:barNum
    node1=deformNode(barConnect(j,1),:);
    node2=deformNode(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color','k');
end




