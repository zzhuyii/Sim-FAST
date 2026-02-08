%% Plot the configuration of the model

function Plot_Shape_Bar_Stress(obj,bar_stress)

% Obtain Information
View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;
assembly=obj.assembly;

%% Set up the graphic space
figure
view(View1,View2); 
set(gca,'DataAspectRatio',[1 1 1])
set(gcf, 'color', 'white');
set(gcf,'position',[obj.x0,obj.y0,obj.width,obj.height])

% The software support two ways to set up the plotting range
A=size(Vsize);
if A(1)==1    
    axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
else
    axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
end

%% Plot the cst element
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
    patch('Faces',f,'Vertices',v,'FaceColor','yellow')
end

%% color code bar stress
barConnect=obj.assembly.bar.node_ij_mat;
barNum=size(barConnect);
barNum=barNum(1);

minSx = min(bar_stress);
maxSx = max(bar_stress);
span  = maxSx - minSx;

for j=1:barNum

    v = bar_stress(j);
    if v > (4/5)*span + minSx
        colorTemp = [1 0 0];        % red
    elseif v > (3/5)*span + minSx
        colorTemp = [1 0.5 0];      % orange (approx RGB)
    elseif v > (2/5)*span + minSx
        colorTemp = [1 1 0];        % yellow
    elseif v > (1/5)*span + minSx
        colorTemp = [0 0.5 0];      % green (darker)
    else
        colorTemp = [0 0 1];        % blue
    end

    node1=assembly.node.coordinates_mat(barConnect(j,1),:);
    node2=assembly.node.coordinates_mat(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color',colorTemp,'LineWidth', 2.0);
end


% Quintile edges from min -> max
edges = minSx + span * (0:5)/5;   % [min, 1/5, 2/5, 3/5, 4/5, max]

% Colors matching your Python code order: red, orange, yellow, green, blue
colors = [ 1    0    0;    % red
           1  0.5    0;    % orange
           1    1    0;    % yellow
           0  0.5    0;    % green
           0    0    1];   % blue

% Bin index pairs for labels (high->low): [4/5,1] [3/5,4/5] [2/5,3/5] [1/5,2/5] [min,1/5]
loIdx = [5 4 3 2 1];
hiIdx = [6 5 4 3 2];

% Create dummy line handles for legend
if ~exist('ax','var') || isempty(ax), ax = gca; end
hold(ax,'on');

h = gobjects(1,5);
labels = cell(1,5);
for k = 1:5
    % Dummy line for legend swatch
    h(k) = plot3(ax, [NaN NaN], [NaN NaN], [NaN NaN], ...
                 'LineWidth', 2, 'Color', colors(k,:));
    % Label in MPa
    labels{k} = sprintf('%.1f to %.1f MPa', ...
        edges(loIdx(k))/1e6, edges(hiIdx(k))/1e6);
end

leg = legend(ax, h, labels, 'Location','northwest');  % upper-left inside
% If you prefer outside like bbox_to_anchor, use:
% leg = legend(ax, h, labels, 'Location','northwestoutside');