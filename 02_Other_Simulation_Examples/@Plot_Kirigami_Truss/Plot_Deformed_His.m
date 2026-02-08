%% plot the deformation history of the simulated results.

function Plot_Deformed_His(obj,Uhis)

% Obtain Information
View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;
assembly=obj.assembly;
undeformedNode=assembly.node.coordinates_mat;
pauseTime=obj.holdTime;
filename=obj.fileName;

%% Set up the graphic space
h=figure;
set(gcf, 'color', 'white');
set(gcf,'position',[obj.x0,obj.y0,obj.width,obj.height])  
B=size(Uhis);
Incre=B(1);

for i=1:Incre
    clf
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    
    A=size(Vsize);
    if A(1)==1    
        axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
    else
        axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
    end
    
    tempU=squeeze(Uhis(i,:,:));
    deformNode=undeformedNode+tempU;

    %% Plot bar
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
    
    %% Plot CST
    cstIJK=obj.assembly.cst.node_ijk_mat;
    panelNum=size(cstIJK);
    panelNum=panelNum(1);
    
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
    
    pause(pauseTime);      

    %% Write to the GIF File 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);    
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 
    end 
end
