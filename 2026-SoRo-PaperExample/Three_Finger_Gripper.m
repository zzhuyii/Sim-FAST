clear all
close all
clc
%% Define Geometry
sectionNum=8;
finger_num=3;

layerL=0.06;
holeL=0.004;
layerH=0.004;

move_dis=0.03; 
% eg. finger move from (0,0,0) to (sqrt(3)*move_dis,move_dis,0), 
% since rotation angle is 60 degree

angle_1=pi/4;    % Rotation by X-axis
angle_2=pi/3;    % Rotation by Z-axis

% film properties
t=0.00011; % 0.11mm for the thickness of film
E=1.3*10^9; % Young's modulus of the film
v=0.2; % Poisson's Ratio of the film

% rotational spring's stiffness
kspr=E*t^3/12;

% target pressure
targetP=30000;

%% Define the Nodal Coordinates
node=Elements_Nodes;

% Finger 1
for i=1:sectionNum
    node.coordinates_mat=[node.coordinates_mat
        -holeL/2 -holeL/2 (i-1)*layerH;
        holeL/2 -holeL/2 (i-1)*layerH;
        holeL/2 holeL/2 (i-1)*layerH;
        -holeL/2 holeL/2 (i-1)*layerH;

        -layerL/4 -layerL/4 (i-0.8)*layerH;
        0 -layerL/4 (i-0.8)*layerH;
        layerL/4 -layerL/4 (i-0.8)*layerH;
        layerL/4 0 (i-0.8)*layerH;
        layerL/4 layerL/4 (i-0.85)*layerH;
        0 layerL/4 (i-0.8)*layerH;
        -layerL/4 layerL/4 (i-0.8)*layerH;
        -layerL/4 0 (i-0.8)*layerH;

        -layerL/2 -layerL/2 (i-0.5)*layerH;
        0 -layerL/2 (i-0.5)*layerH;
        layerL/2 -layerL/2 (i-0.5)*layerH;
        layerL/2 0 (i-0.5)*layerH;
        layerL/2 layerL/2 (i-0.5)*layerH;
        0 layerL/2 (i-0.5)*layerH;
        -layerL/2 layerL/2 (i-0.5)*layerH;
        -layerL/2 0 (i-0.5)*layerH;

        -layerL/4 -layerL/4 (i-0.2)*layerH;
        0 -layerL/4 (i-0.2)*layerH;
        layerL/4 -layerL/4 (i-0.2)*layerH;
        layerL/4 0 (i-0.2)*layerH;
        layerL/4 layerL/4 (i-0.2)*layerH;
        0 layerL/4 (i-0.2)*layerH;
        -layerL/4 layerL/4 (i-0.2)*layerH;
        -layerL/4 0 (i-0.2)*layerH;
        ];
end

node.coordinates_mat=[node.coordinates_mat
    -holeL/2 -holeL/2 sectionNum*layerH;
    holeL/2 -holeL/2 sectionNum*layerH;
    holeL/2 holeL/2 sectionNum*layerH;
    -holeL/2 holeL/2 sectionNum*layerH];

% Number for nodes in first finger 
finger_1_node_num=size(node.coordinates_mat,1); 

% Rotate the finger to the appropriate location
% first rotation by X-axis
RotMat=[1 0 0;
        0 cos(-angle_1) sin(-angle_1);
        0 -sin(-angle_1) cos(-angle_1)]; 
node.coordinates_mat=node.coordinates_mat*RotMat;

% second rotation by Z-axis
RotMat=[cos(-angle_2)   sin(-angle_2) 0;
        -sin(-angle_2)  cos(-angle_2) 0;
        0 0 1];
node.coordinates_mat=node.coordinates_mat*RotMat; 

% Move the finger to the appropriate location
node.coordinates_mat(:,1)=node.coordinates_mat(:,1)+sqrt(3)*move_dis;
node.coordinates_mat(:,2)=node.coordinates_mat(:,2)+move_dis;

% Finger 2
for i=1:sectionNum
    node.coordinates_mat=[node.coordinates_mat
        -holeL/2 -holeL/2 (i-1)*layerH;
        holeL/2 -holeL/2 (i-1)*layerH;
        holeL/2 holeL/2 (i-1)*layerH;
        -holeL/2 holeL/2 (i-1)*layerH;

        -layerL/4 -layerL/4 (i-0.8)*layerH;
        0 -layerL/4 (i-0.8)*layerH;
        layerL/4 -layerL/4 (i-0.8)*layerH;
        layerL/4 0 (i-0.8)*layerH;
        layerL/4 layerL/4 (i-0.85)*layerH;
        0 layerL/4 (i-0.8)*layerH;
        -layerL/4 layerL/4 (i-0.8)*layerH;
        -layerL/4 0 (i-0.8)*layerH;

        -layerL/2 -layerL/2 (i-0.5)*layerH;
        0 -layerL/2 (i-0.5)*layerH;
        layerL/2 -layerL/2 (i-0.5)*layerH;
        layerL/2 0 (i-0.5)*layerH;
        layerL/2 layerL/2 (i-0.5)*layerH;
        0 layerL/2 (i-0.5)*layerH;
        -layerL/2 layerL/2 (i-0.5)*layerH;
        -layerL/2 0 (i-0.5)*layerH;

        -layerL/4 -layerL/4 (i-0.2)*layerH;
        0 -layerL/4 (i-0.2)*layerH;
        layerL/4 -layerL/4 (i-0.2)*layerH;
        layerL/4 0 (i-0.2)*layerH;
        layerL/4 layerL/4 (i-0.2)*layerH;
        0 layerL/4 (i-0.2)*layerH;
        -layerL/4 layerL/4 (i-0.2)*layerH;
        -layerL/4 0 (i-0.2)*layerH;
        ];
end

node.coordinates_mat=[node.coordinates_mat
    -holeL/2 -holeL/2 sectionNum*layerH;
    holeL/2 -holeL/2 sectionNum*layerH;
    holeL/2 holeL/2 sectionNum*layerH;
    -holeL/2 holeL/2 sectionNum*layerH];

% Node number of second finger
finger_2_node_num=size(node.coordinates_mat,1)-finger_1_node_num; 

% Rotate the finger to the appropriate location
% first rotation by X-axis
RotMat=[1 0 0;
        0 cos(-angle_1) sin(-angle_1);
        0 -sin(-angle_1) cos(-angle_1)]; 
node.coordinates_mat(finger_1_node_num+1:end,:)=node.coordinates_mat(finger_1_node_num+1:end,:)*RotMat;

% Rotation by Z-axis
RotMat=[cos(angle_2)   sin(angle_2) 0;
        -sin(angle_2)  cos(angle_2) 0;
        0 0 1];
node.coordinates_mat(finger_1_node_num+1:end,:)=node.coordinates_mat(finger_1_node_num+1:end,:)*RotMat; 

% Move the finger to the appropriate location
node.coordinates_mat(finger_1_node_num+1:end,1)=node.coordinates_mat(finger_1_node_num+1:end,1)-sqrt(3)*move_dis;
node.coordinates_mat(finger_1_node_num+1:end,2)=node.coordinates_mat(finger_1_node_num+1:end,2)+move_dis;

% Finger 3
for i=1:sectionNum
    node.coordinates_mat=[node.coordinates_mat
        -holeL/2 -holeL/2 (i-1)*layerH;
        holeL/2 -holeL/2 (i-1)*layerH;
        holeL/2 holeL/2 (i-1)*layerH;
        -holeL/2 holeL/2 (i-1)*layerH;

        -layerL/4 -layerL/4 (i-0.8)*layerH;
        0 -layerL/4 (i-0.8)*layerH;
        layerL/4 -layerL/4 (i-0.8)*layerH;
        layerL/4 0 (i-0.8)*layerH;
        layerL/4 layerL/4 (i-0.85)*layerH;
        0 layerL/4 (i-0.8)*layerH;
        -layerL/4 layerL/4 (i-0.8)*layerH;
        -layerL/4 0 (i-0.8)*layerH;

        -layerL/2 -layerL/2 (i-0.5)*layerH;
        0 -layerL/2 (i-0.5)*layerH;
        layerL/2 -layerL/2 (i-0.5)*layerH;
        layerL/2 0 (i-0.5)*layerH;
        layerL/2 layerL/2 (i-0.5)*layerH;
        0 layerL/2 (i-0.5)*layerH;
        -layerL/2 layerL/2 (i-0.5)*layerH;
        -layerL/2 0 (i-0.5)*layerH;

        -layerL/4 -layerL/4 (i-0.2)*layerH;
        0 -layerL/4 (i-0.2)*layerH;
        layerL/4 -layerL/4 (i-0.2)*layerH;
        layerL/4 0 (i-0.2)*layerH;
        layerL/4 layerL/4 (i-0.2)*layerH;
        0 layerL/4 (i-0.2)*layerH;
        -layerL/4 layerL/4 (i-0.2)*layerH;
        -layerL/4 0 (i-0.2)*layerH;
        ];
end

node.coordinates_mat=[node.coordinates_mat
    -holeL/2 -holeL/2 sectionNum*layerH;
    holeL/2 -holeL/2 sectionNum*layerH;
    holeL/2 holeL/2 sectionNum*layerH;
    -holeL/2 holeL/2 sectionNum*layerH];

% the node number for the third finger 
finger_3_node_num=size(node.coordinates_mat,1)-finger_1_node_num-finger_2_node_num; 

% Rotate the finger to the appropriate location
% first rotation by X-axis
RotMat=[1 0 0;
        0 cos(angle_1) sin(angle_1);
        0 -sin(angle_1) cos(angle_1)]; 
node.coordinates_mat(finger_1_node_num+finger_2_node_num+1:end,:)=node.coordinates_mat(finger_1_node_num+finger_2_node_num+1:end,:)*RotMat;

% Move the finger to the right location
node.coordinates_mat(finger_1_node_num+finger_2_node_num+1:end,2)=node.coordinates_mat(finger_1_node_num+finger_2_node_num+1:end,2)-2*move_dis;

%% Define assembly
assembly=Assembly_Membrane;
assembly.node=node;

cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings_4N;

assembly.cst=cst;
assembly.rotSpr=rotSpr;

%% Define Plotting Functions
plots=Plot_Membrane;
plots.assembly=assembly;
plots.displayRange=[-3; 3; -3; 3; -3; sectionNum]*layerL;

% plot the nodal location for inspection
plots.Plot_Shape_NodeNumber;

%% Set up the Triangle Elements
tri_ijk=[];
tri_direction=[];
each_finger_node=finger_1_node_num;

for j=1:finger_num
    for i=1:sectionNum  
        tri_ijk=[tri_ijk;
            (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node;    
            (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node;   
            (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node;   
            (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node;
    
            (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node;    
            (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node (i-1)*28+7+(j-1)*each_finger_node;   
            (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node;   
            (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node (i-1)*28+9+(j-1)*each_finger_node;
            (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node; 
            (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node (i-1)*28+11+(j-1)*each_finger_node; 
            (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node;
            (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node (i-1)*28+5+(j-1)*each_finger_node; 
    
            (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+13+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;    
            (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;   
            (i-1)*28+6+(j-1)*each_finger_node (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;   
            (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node (i-1)*28+15+(j-1)*each_finger_node;
            (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+15+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node; 
            (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node; 
            (i-1)*28+8+(j-1)*each_finger_node (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node;
            (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node (i-1)*28+17+(j-1)*each_finger_node; 
    
            (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+17+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;    
            (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;   
            (i-1)*28+10+(j-1)*each_finger_node (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;   
            (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node (i-1)*28+19+(j-1)*each_finger_node;
            (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+19+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node; 
            (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node; 
            (i-1)*28+12+(j-1)*each_finger_node (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node;
            (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node (i-1)*28+13+(j-1)*each_finger_node; 
            ];   
    
        % The bottom film is facing down
        tri_direction=[tri_direction;
            -ones(28,1)];
    
        tri_ijk=[tri_ijk;
            (i)*28+1+(j-1)*each_finger_node (i)*28+2+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node;    
            (i)*28+2+(j-1)*each_finger_node (i)*28+3+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node;   
            (i)*28+3+(j-1)*each_finger_node (i)*28+4+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node;   
            (i)*28+4+(j-1)*each_finger_node (i)*28+1+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node;
    
            (i)*28+1+(j-1)*each_finger_node (i-1)*28+5+16+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node;    
            (i)*28+2+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node (i-1)*28+7+16+(j-1)*each_finger_node;   
            (i)*28+2+(j-1)*each_finger_node (i-1)*28+7+16+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node;   
            (i)*28+3+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node (i-1)*28+9+16+(j-1)*each_finger_node;
            (i)*28+3+(j-1)*each_finger_node (i-1)*28+9+16+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node; 
            (i)*28+4+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node (i-1)*28+11+16+(j-1)*each_finger_node; 
            (i)*28+4+(j-1)*each_finger_node (i-1)*28+11+16+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node;
            (i)*28+1+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node (i-1)*28+5+16+(j-1)*each_finger_node; 
    
            (i-1)*28+5+16+(j-1)*each_finger_node (i-1)*28+13+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;    
            (i-1)*28+5+16+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;   
            (i-1)*28+6+16+(j-1)*each_finger_node (i-1)*28+7+16+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;   
            (i-1)*28+7+16+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node (i-1)*28+15+(j-1)*each_finger_node;
            (i-1)*28+7+16+(j-1)*each_finger_node (i-1)*28+15+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node; 
            (i-1)*28+7+16+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node; 
            (i-1)*28+8+16+(j-1)*each_finger_node (i-1)*28+9+16+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node;
            (i-1)*28+9+16+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node (i-1)*28+17+(j-1)*each_finger_node; 
    
            (i-1)*28+9+16+(j-1)*each_finger_node (i-1)*28+17+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;    
            (i-1)*28+9+16+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;   
            (i-1)*28+10+16+(j-1)*each_finger_node (i-1)*28+11+16+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;   
            (i-1)*28+11+16+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node (i-1)*28+19+(j-1)*each_finger_node;
            (i-1)*28+11+16+(j-1)*each_finger_node (i-1)*28+19+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node; 
            (i-1)*28+11+16+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node; 
            (i-1)*28+12+16+(j-1)*each_finger_node (i-1)*28+5+16+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node;
            (i-1)*28+5+16+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node (i-1)*28+13+(j-1)*each_finger_node; 
            ];   
    
        % The top film is facing up
        tri_direction=[tri_direction;
            ones(28,1)];
    end
end

% number of triangles
triNum=size(tri_ijk,1);

% Organize the trinagle sequence so that 
% the normal vector points to the out side
for i=1:triNum
    x1=node.coordinates_mat(tri_ijk(i,1),:);
    x2=node.coordinates_mat(tri_ijk(i,2),:);
    x3=node.coordinates_mat(tri_ijk(i,3),:);
    v1=x2-x1;
    v2=x3-x1;
    normal=cross(v1,v2);
    if sign(normal(3))==tri_direction(i)
    else
        temp=tri_ijk(i,2);
        tri_ijk(i,2)=tri_ijk(i,3);
        tri_ijk(i,3)=temp;
    end
end

% Define Triangle Elements
cst.node_ijk_mat=tri_ijk;

% Define Young's modulus, thickness, Poisson's ratio
cst.E_vec=E*ones(triNum,1);
cst.t_vec=t*ones(triNum,1);
cst.v_vec=v*ones(triNum,1);

% Plot the CST elements for inspection
plots.Plot_Shape_CSTNumber;

%% Define Rotational Spring
for j=1:finger_num
    for i=1:sectionNum
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;    
            (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;
            (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;
            (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node;
            (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node;
            (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;
            (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;
            (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node;
            (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node;
    
            (i)*28+1+(j-1)*each_finger_node (i-1)*28+5+16+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;
            (i)*28+2+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node (i-1)*28+7+16+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node;
            (i)*28+2+(j-1)*each_finger_node (i-1)*28+7+16+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node;
            (i)*28+3+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node (i-1)*28+9+16+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node;
            (i)*28+3+(j-1)*each_finger_node (i-1)*28+9+16+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;
            (i)*28+4+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node (i-1)*28+11+16+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node;
            (i)*28+4+(j-1)*each_finger_node (i-1)*28+11+16+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node;
            (i)*28+1+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node (i-1)*28+5+16+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node;
    
            (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node (i-1)*28+7+(j-1)*each_finger_node;
            (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node (i-1)*28+9+(j-1)*each_finger_node;
            (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node (i-1)*28+11+(j-1)*each_finger_node;
            (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node (i-1)*28+5+(j-1)*each_finger_node;
    
            (i-1)*28+5+16+(j-1)*each_finger_node (i-1)*28+14+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node (i-1)*28+7+16+(j-1)*each_finger_node;
            (i-1)*28+7+16+(j-1)*each_finger_node (i-1)*28+16+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node (i-1)*28+9+16+(j-1)*each_finger_node;
            (i-1)*28+9+16+(j-1)*each_finger_node (i-1)*28+18+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node (i-1)*28+11+16+(j-1)*each_finger_node;
            (i-1)*28+11+16+(j-1)*each_finger_node (i-1)*28+20+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node (i-1)*28+5+16+(j-1)*each_finger_node;
    
            (i-1)*28+5+(j-1)*each_finger_node (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node (i-1)*28+2+(j-1)*each_finger_node;
            (i-1)*28+7+(j-1)*each_finger_node (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node (i-1)*28+3+(j-1)*each_finger_node;
            (i-1)*28+9+(j-1)*each_finger_node (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node (i-1)*28+4+(j-1)*each_finger_node;
            (i-1)*28+11+(j-1)*each_finger_node (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node (i-1)*28+1+(j-1)*each_finger_node;
    
            (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+6+(j-1)*each_finger_node (i-1)*28+7+(j-1)*each_finger_node;
            (i-1)*28+2+(j-1)*each_finger_node (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+8+(j-1)*each_finger_node (i-1)*28+9+(j-1)*each_finger_node;
            (i-1)*28+3+(j-1)*each_finger_node (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+10+(j-1)*each_finger_node (i-1)*28+11+(j-1)*each_finger_node;
            (i-1)*28+4+(j-1)*each_finger_node (i-1)*28+1+(j-1)*each_finger_node (i-1)*28+12+(j-1)*each_finger_node (i-1)*28+5+(j-1)*each_finger_node;
    
            (i-1)*28+5+16+(j-1)*each_finger_node (i)*28+1+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node (i)*28+2+(j-1)*each_finger_node;
            (i-1)*28+7+16+(j-1)*each_finger_node (i)*28+2+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node (i)*28+3+(j-1)*each_finger_node;
            (i-1)*28+9+16+(j-1)*each_finger_node (i)*28+3+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node (i)*28+4+(j-1)*each_finger_node;
            (i-1)*28+11+16+(j-1)*each_finger_node (i)*28+4+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node (i)*28+1+(j-1)*each_finger_node;
    
            (i)*28+1+(j-1)*each_finger_node (i)*28+2+(j-1)*each_finger_node (i-1)*28+6+16+(j-1)*each_finger_node (i-1)*28+7+16+(j-1)*each_finger_node;
            (i)*28+2+(j-1)*each_finger_node (i)*28+3+(j-1)*each_finger_node (i-1)*28+8+16+(j-1)*each_finger_node (i-1)*28+9+16+(j-1)*each_finger_node;
            (i)*28+3+(j-1)*each_finger_node (i)*28+4+(j-1)*each_finger_node (i-1)*28+10+16+(j-1)*each_finger_node (i-1)*28+11+16+(j-1)*each_finger_node;
            (i)*28+4+(j-1)*each_finger_node (i)*28+1+(j-1)*each_finger_node (i-1)*28+12+16+(j-1)*each_finger_node (i-1)*28+5+16+(j-1)*each_finger_node;
            ];
        if i<sectionNum
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*28+22+(j-1)*each_finger_node (i-1)*28+29+(j-1)*each_finger_node (i-1)*28+30+(j-1)*each_finger_node (i-1)*28+34+(j-1)*each_finger_node;
                (i-1)*28+24+(j-1)*each_finger_node (i-1)*28+30+(j-1)*each_finger_node (i-1)*28+31+(j-1)*each_finger_node (i-1)*28+36+(j-1)*each_finger_node;
                (i-1)*28+26+(j-1)*each_finger_node (i-1)*28+31+(j-1)*each_finger_node (i-1)*28+32+(j-1)*each_finger_node (i-1)*28+38+(j-1)*each_finger_node;
                (i-1)*28+28+(j-1)*each_finger_node (i-1)*28+32+(j-1)*each_finger_node (i-1)*28+29+(j-1)*each_finger_node (i-1)*28+40+(j-1)*each_finger_node;
                ];
        end
    end
end

rotSpr.theta1=0.02*pi;
rotSpr.theta2=(2-0.02)*pi;

% Find the bending stiffness
sprNum=size(rotSpr.node_ijkl_mat,1);
rotSpr.rot_spr_K_vec=kspr*ones(sprNum,1);

% Plot the rotational springs for inspection
plots.Plot_Shape_SprNumber;

%% Initialize the assembly
assembly.Initialize_Assembly;

%% Set up solver for inflation
nr=Solver_NR_Loading;
nr.assembly=assembly;

for i=1:finger_num
    nr.supp=[nr.supp;
        13+(i-1)*each_finger_node 1 1 1;    
        15+(i-1)*each_finger_node 1 1 1;
        17+(i-1)*each_finger_node 1 1 1;
        19+(i-1)*each_finger_node 1 1 1;
    ];
end

step=150;
Uhis=[];

for k=1:step

    % To help with convergence, we further make the first few steps 
    % to have smaller step lengths 
    if k<3
        subStep=30;
        for q=1:subStep
            % Update the pressure forces for nonlinearity
            xcurrent=node.coordinates_mat+node.current_U_mat;    
            nodeNum=size(node.coordinates_mat,1);
            fmat=SolvePressureForce(xcurrent,(k-1+q/subStep)/step*targetP,tri_ijk);
            fpressure=[(1:nodeNum)',fmat];
        
            nr.load=fpressure;            
            nr.increStep=1;
            nr.iterMax=20;
            nr.tol=10^-3;
        
            Utemp=nr.Solve();
        end
        Uhis(k,:,:)=squeeze(nr.Solve());
    end

    % Update the pressure forces for nonlinearity
    xcurrent=node.coordinates_mat+node.current_U_mat;
    nodeNum=size(node.coordinates_mat,1);
    fmat=SolvePressureForce(xcurrent,k/step*targetP,tri_ijk);
    fpressure=[(1:nodeNum)',fmat];

    nr.load=fpressure;    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=10^-6;

    Uhis(k,:,:)=squeeze(nr.Solve());

end

% Plot the innitial configuration
figure;
plots.Plot_DeformedShape(squeeze(Uhis(1,:,:)))

% Plot the inflated configuration
figure;
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
UhisExtend=Uhis;


%% Set up the solver for finger bending
for i=1:finger_num
nr.supp=[nr.supp;
    13+(i-1)*each_finger_node 1 1 1;    
    15+(i-1)*each_finger_node 1 1 1;
    17+(i-1)*each_finger_node 1 1 1;
    19+(i-1)*each_finger_node 1 1 1;
];
end

step=100;
targetF=3;
frictionLoss=0.35;
Uhis=[];

for k=1:step

    % Update the pressure forces for nonlinearity
    xcurrent=node.coordinates_mat+node.current_U_mat;
    nodeNum=size(node.coordinates_mat,1);
    fmat=SolvePressureForce(xcurrent,targetP,tri_ijk);
    fpressure=[(1:nodeNum)',fmat];

    for i=1:sectionNum-1
    
        n1=xcurrent((i)*28+13,:);
        n2=xcurrent((i-1)*28+13,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);
   
        fpressure((i)*28+13,2:4)=fpressure((i)*28+13,2:4)+v1*(targetF-frictionLoss*(i-1))*k/step;
        fpressure((i-1)*28+13,2:4)=fpressure((i-1)*28+13,2:4)+v2*(targetF-frictionLoss*(i-1))*k/step;

        n1=xcurrent((i)*28+15,:);
        n2=xcurrent((i-1)*28+15,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*28+15,2:4)=fpressure((i)*28+15,2:4)+v1*(targetF-frictionLoss*(i-1))*k/step;
        fpressure((i-1)*28+15,2:4)=fpressure((i-1)*28+15,2:4)+v2*(targetF-frictionLoss*(i-1))*k/step; 

        n1=xcurrent((i)*28+13+each_finger_node,:);
        n2=xcurrent((i-1)*28+13+each_finger_node,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*28+13+each_finger_node,2:4)=fpressure((i)*28+13+each_finger_node,2:4)+v1*(targetF-frictionLoss*(i-1))*k/step;
        fpressure((i-1)*28+13+each_finger_node,2:4)=fpressure((i-1)*28+13+each_finger_node,2:4)+v2*(targetF-frictionLoss*(i-1))*k/step;

        n1=xcurrent((i)*28+15+each_finger_node,:);
        n2=xcurrent((i-1)*28+15+each_finger_node,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*28+15+each_finger_node,2:4)=fpressure((i)*28+15+each_finger_node,2:4)+v1*(targetF-frictionLoss*(i-1))*k/step;
        fpressure((i-1)*28+15+each_finger_node,2:4)=fpressure((i-1)*28+15+each_finger_node,2:4)+v2*(targetF-frictionLoss*(i-1))*k/step; % finger 2 bending

        n1=xcurrent((i)*28+17+2*each_finger_node,:);
        n2=xcurrent((i-1)*28+17+2*each_finger_node,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*28+17+2*each_finger_node,2:4)=fpressure((i)*28+17+2*each_finger_node,2:4)+v1*(targetF-frictionLoss*(i-1))*k/step;
        fpressure((i-1)*28+17+2*each_finger_node,2:4)=fpressure((i-1)*28+17+2*each_finger_node,2:4)+v2*(targetF-frictionLoss*(i-1))*k/step;

        n1=xcurrent((i)*28+19+2*each_finger_node,:);
        n2=xcurrent((i-1)*28+19+2*each_finger_node,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*28+19+2*each_finger_node,2:4)=fpressure((i)*28+19+2*each_finger_node,2:4)+v1*(targetF-frictionLoss*(i-1))*k/step;
        fpressure((i-1)*28+19+2*each_finger_node,2:4)=fpressure((i-1)*28+19+2*each_finger_node,2:4)+v2*(targetF-frictionLoss*(i-1))*k/step; % finger 3 bending

    end

    nr.load=fpressure;    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=10^-6;

    Uhis(k,:,:)=squeeze(nr.Solve());

end

% Plot the bending configuration
figure;
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))

% Get the full animation
plots.fileName='Three_Finger_Gripper.gif';
Uhis=cat(1,UhisExtend,Uhis);
plots.Plot_DeformedHis(Uhis(1:10:end,:,:))

%% Solve the applied force due to pressure
function Fmat=SolvePressureForce(NodeCord,P,tri_ijk)
    nodeNum=size(NodeCord,1);
    Fmat=zeros(nodeNum,3);
    triNum=size(tri_ijk,1);

    for i=1:triNum
        x1=NodeCord(tri_ijk(i,1),:);
        x2=NodeCord(tri_ijk(i,2),:);
        x3=NodeCord(tri_ijk(i,3),:);

        v1=x2-x1;
        v2=x3-x1;

        normal=cross(v1,v2);
        A=norm(normal)/2;

        f=A*P/3;
        normal=normal/norm(normal);

        Fmat(tri_ijk(i,1),:)=Fmat(tri_ijk(i,1),:)+f*normal;
        Fmat(tri_ijk(i,2),:)=Fmat(tri_ijk(i,2),:)+f*normal;
        Fmat(tri_ijk(i,3),:)=Fmat(tri_ijk(i,3),:)+f*normal;
    end   
end