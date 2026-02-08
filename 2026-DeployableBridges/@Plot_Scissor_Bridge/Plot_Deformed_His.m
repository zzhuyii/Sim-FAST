%% Plot the deformation history of the simulated results.
function Plot_Deformed_His(obj, Uhis)

% Obtain Information
View1 = obj.viewAngle1;
View2 = obj.viewAngle2;
Vsize = obj.displayRange;
Vratio = obj.displayRangeRatio;
assembly = obj.assembly;
undeformedNode = assembly.node.coordinates_mat;  
pauseTime = obj.holdTime;
filename = obj.fileName;

%% Set up the graphic space
h = figure;
set(gcf, 'color', 'white');
set(gcf, 'position', [obj.x0, obj.y0, obj.width, obj.height]);
Incre = size(Uhis,1);

for i = 1:Incre
    clf
    view(View1, View2);
    set(gca, 'DataAspectRatio', [1 1 1])
    
    if numel(Vsize) == 2
        axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
    else
        axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
    end    
    tempU = squeeze(Uhis(i,:,:));
    deformNode = undeformedNode + tempU;

    %% Plot the bar element
    barConnect = obj.assembly.bar.node_ij_mat;
    barNum = size(barConnect,1);
    for j = 1:barNum
        node1 = deformNode(barConnect(j,1), :);
        node2 = deformNode(barConnect(j,2), :);
        line([node1(1), node2(1)], ...
             [node1(2), node2(2)], ...
             [node1(3), node2(3)], 'Color', 'k');
    end


    %% Plot the active bar
    actBarConnect=assembly.actBar.node_ij_mat;
    actBarNum=size(actBarConnect);
    actBarNum=actBarNum(1);

    for j=1:actBarNum
        node1=deformNode(actBarConnect(j,1),:);
        node2=deformNode(actBarConnect(j,2),:);
        line([node1(1),node2(1)],...
             [node1(2),node2(2)],...
             [node1(3),node2(3)],'Color','b','Linewidth',3);
    end
    

    %% Plot the CST element
    cstIJK = obj.assembly.cst.node_ijk_mat;
    panelNum = size(cstIJK,1);
    for k = 1:panelNum
        nodeNumVec = cstIJK(k, :);
        f = 1:length(nodeNumVec);
        v = zeros(length(nodeNumVec), 3);
        for j = 1:length(nodeNumVec)
            v(j, :) = deformNode(nodeNumVec(j), :);
        end
        patch('Faces', f, 'Vertices', v, 'FaceColor', 'yellow');
    end
    
    %% Generate the GIF figure
    pause(pauseTime);    
    frame = getframe(h);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', pauseTime);
    end
end
