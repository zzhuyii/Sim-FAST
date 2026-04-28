clear all
close all
clc

%% Define Geometry

% Section amount
sectionNum=10;

% Bottom layer hexagon side length & Shrinking rate from bottom to top
Bottom_Hex_side_length=0.05;
side_shrink_rate=0.95;

% Hexagon side length on each layer
Hex_side_length=[];
for i=1:sectionNum
Hex_side_length=[Hex_side_length
         Bottom_Hex_side_length*(side_shrink_rate^(i-1));];
end

% Internal tunnel size (Length of hole hexagon side)
Hole_hex_side_length=0.01;

% Layer height initially
layerH=0.002;

% Plastic films properties
t=0.00011; % Thickness 0.11mm
E=1.3*10^9; % Young's modulus
v=0.2; % Possion ratio

% Stiffness of rotational spring
kspr=10*E*t^3/12;

% nodes initial coordinates of hexagon structure
node=Elements_Nodes;
for i=1:sectionNum
    node.coordinates_mat=[node.coordinates_mat
        Hole_hex_side_length 0 (i-1)*layerH;
        (1/2)*Hole_hex_side_length (sqrt(3)/2)*Hole_hex_side_length (i-1)*layerH;
        -(1/2)*Hole_hex_side_length (sqrt(3)/2)*Hole_hex_side_length (i-1)*layerH;
        -Hole_hex_side_length 0 (i-1)*layerH;
        -(1/2)*Hole_hex_side_length -(sqrt(3)/2)*Hole_hex_side_length (i-1)*layerH;
        (1/2)*Hole_hex_side_length -(sqrt(3)/2)*Hole_hex_side_length (i-1)*layerH; % inside nodes (hole)

        (Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length 0 (i-1)*layerH+(1/4)*layerH;
        (3/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (sqrt(3)/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        (1/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        0 (sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        -(1/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        -(3/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (sqrt(3)/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        -((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) 0 (i-1)*layerH+(1/4)*layerH;
        -(3/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) -(sqrt(3)/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        -(1/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) -(sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        0 -(sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        (1/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) -(sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH;
        (3/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) -(sqrt(3)/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(1/4)*layerH; % middle point 1st layer
        
        Hex_side_length(i) 0 (i-1)*layerH+(1/2)*layerH;
        (3/4)*Hex_side_length(i) (sqrt(3)/4)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        (1/2)*Hex_side_length(i) (sqrt(3)/2)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        0 (sqrt(3)/2)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        -(1/2)*Hex_side_length(i) (sqrt(3)/2)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        -(3/4)*Hex_side_length(i) (sqrt(3)/4)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        -Hex_side_length(i) 0 (i-1)*layerH+(1/2)*layerH;
        -(3/4)*Hex_side_length(i) -(sqrt(3)/4)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        -(1/2)*Hex_side_length(i) -(sqrt(3)/2)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        0 -(sqrt(3)/2)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        (1/2)*Hex_side_length(i) -(sqrt(3)/2)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;
        (3/4)*Hex_side_length(i) -(sqrt(3)/4)*Hex_side_length(i) (i-1)*layerH+(1/2)*layerH;% outside nodes (Big hexagon)

        (Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length 0 (i-1)*layerH+(3/4)*layerH;
        (3/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (sqrt(3)/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        (1/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        0 (sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        -(1/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        -(3/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (sqrt(3)/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        -((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) 0 (i-1)*layerH+(3/4)*layerH;
        -(3/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) -(sqrt(3)/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        -(1/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) -(sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        0 -(sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        (1/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) -(sqrt(3)/2)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH;
        (3/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) -(sqrt(3)/4)*((Hex_side_length(i)-Hole_hex_side_length)/2+Hole_hex_side_length) (i-1)*layerH+(3/4)*layerH; % middle point 2st layer

        ];
end
node.coordinates_mat=[node.coordinates_mat
        Hole_hex_side_length 0 sectionNum*layerH;
        (1/2)*Hole_hex_side_length (sqrt(3)/2)*Hole_hex_side_length sectionNum*layerH;
        -(1/2)*Hole_hex_side_length (sqrt(3)/2)*Hole_hex_side_length sectionNum*layerH;
        -Hole_hex_side_length 0 sectionNum*layerH;
        -(1/2)*Hole_hex_side_length -(sqrt(3)/2)*Hole_hex_side_length sectionNum*layerH;
        (1/2)*Hole_hex_side_length -(sqrt(3)/2)*Hole_hex_side_length sectionNum*layerH;]; % Top cap inside nodes

% Rotate about z axis
angle_1=pi/6;    % Rotation by X-axis
RotMat=[cos(-angle_1) sin(-angle_1) 0;
        -sin(-angle_1) cos(-angle_1) 0;
        0 0 1;]; % Rotation by X-axis
node.coordinates_mat=node.coordinates_mat*RotMat;

% Relocate whole hexagon structure
% Also reflect the robot to face bottom
moving_dis = 0.3;  % move 0.5 meters in +X direction
node.coordinates_mat(:,1) = node.coordinates_mat(:,1) + moving_dis;
moving_dis = 0.7;  % move 0.5 meters in +Z direction
node.coordinates_mat(:,3) = -node.coordinates_mat(:,3) + moving_dis;

%% Define assembly
assembly=Assembly_Membrane;
assembly.node=node;

cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings_4N;
actBar=CD_Elements_Bars;

% assembly.actBar=actBar;
assembly.cst=cst;
assembly.rotSpr=rotSpr;

%% Define Plotting Functions
plots=Plot_Membrane;
plots.assembly=assembly;
plots.displayRange=[0.2; 0.9; -0.2; 0.2; -0.02; 0.42];

plots.Plot_Shape_NodeNumber;

%% Set up the triangle information for pressure
tri_ijk=[];
tri_direction=[];
for i=1:sectionNum  
    tri_ijk=[tri_ijk;
        (i-1)*42+1 (i-1)*42+2 (i-1)*42+8;
        (i-1)*42+2 (i-1)*42+3 (i-1)*42+10;
        (i-1)*42+3 (i-1)*42+4 (i-1)*42+12;
        (i-1)*42+4 (i-1)*42+5 (i-1)*42+14;
        (i-1)*42+5 (i-1)*42+6 (i-1)*42+16;
        (i-1)*42+6 (i-1)*42+1 (i-1)*42+18;

        (i-1)*42+1 (i-1)*42+7 (i-1)*42+8;
        (i-1)*42+2 (i-1)*42+8 (i-1)*42+9;
        (i-1)*42+2 (i-1)*42+9 (i-1)*42+10;
        (i-1)*42+3 (i-1)*42+10 (i-1)*42+11;
        (i-1)*42+3 (i-1)*42+11 (i-1)*42+12;
        (i-1)*42+4 (i-1)*42+12 (i-1)*42+13;
        (i-1)*42+4 (i-1)*42+13 (i-1)*42+14;
        (i-1)*42+5 (i-1)*42+14 (i-1)*42+15;
        (i-1)*42+5 (i-1)*42+15 (i-1)*42+16;
        (i-1)*42+6 (i-1)*42+16 (i-1)*42+17;
        (i-1)*42+6 (i-1)*42+17 (i-1)*42+18;
        (i-1)*42+1 (i-1)*42+18 (i-1)*42+7;

        (i-1)*42+7 (i-1)*42+8 (i-1)*42+20;
        (i-1)*42+8 (i-1)*42+9 (i-1)*42+20;
        (i-1)*42+9 (i-1)*42+10 (i-1)*42+22;
        (i-1)*42+10 (i-1)*42+11 (i-1)*42+22;
        (i-1)*42+11 (i-1)*42+12 (i-1)*42+24;
        (i-1)*42+12 (i-1)*42+13 (i-1)*42+24;
        (i-1)*42+13 (i-1)*42+14 (i-1)*42+26;
        (i-1)*42+14 (i-1)*42+15 (i-1)*42+26;
        (i-1)*42+15 (i-1)*42+16 (i-1)*42+28;
        (i-1)*42+16 (i-1)*42+17 (i-1)*42+28;
        (i-1)*42+17 (i-1)*42+18 (i-1)*42+30;
        (i-1)*42+18 (i-1)*42+7 (i-1)*42+30;

        (i-1)*42+7 (i-1)*42+19 (i-1)*42+20;
        (i-1)*42+9 (i-1)*42+20 (i-1)*42+21;
        (i-1)*42+9 (i-1)*42+21 (i-1)*42+22;
        (i-1)*42+11 (i-1)*42+22 (i-1)*42+23;
        (i-1)*42+11 (i-1)*42+23 (i-1)*42+24;
        (i-1)*42+13 (i-1)*42+24 (i-1)*42+25;
        (i-1)*42+13 (i-1)*42+25 (i-1)*42+26;
        (i-1)*42+15 (i-1)*42+26 (i-1)*42+27;
        (i-1)*42+15 (i-1)*42+27 (i-1)*42+28;
        (i-1)*42+17 (i-1)*42+28 (i-1)*42+29;
        (i-1)*42+17 (i-1)*42+29 (i-1)*42+30;
        (i-1)*42+7 (i-1)*42+30 (i-1)*42+19;

        ];
    tri_direction=[tri_direction;
        -ones(42,1)];

    tri_ijk=[tri_ijk;

        (i-1)*42+31 (i-1)*42+19 (i-1)*42+20;
        (i-1)*42+33 (i-1)*42+20 (i-1)*42+21;
        (i-1)*42+33 (i-1)*42+21 (i-1)*42+22;
        (i-1)*42+35 (i-1)*42+22 (i-1)*42+23;
        (i-1)*42+35 (i-1)*42+23 (i-1)*42+24;
        (i-1)*42+37 (i-1)*42+24 (i-1)*42+25;
        (i-1)*42+37 (i-1)*42+25 (i-1)*42+26;
        (i-1)*42+39 (i-1)*42+26 (i-1)*42+27;
        (i-1)*42+39 (i-1)*42+27 (i-1)*42+28;
        (i-1)*42+41 (i-1)*42+28 (i-1)*42+29;
        (i-1)*42+41 (i-1)*42+29 (i-1)*42+30;
        (i-1)*42+31 (i-1)*42+30 (i-1)*42+19;

        (i-1)*42+31 (i-1)*42+32 (i-1)*42+20;
        (i-1)*42+32 (i-1)*42+33 (i-1)*42+20;
        (i-1)*42+33 (i-1)*42+34 (i-1)*42+22;
        (i-1)*42+34 (i-1)*42+35 (i-1)*42+22;
        (i-1)*42+35 (i-1)*42+36 (i-1)*42+24;
        (i-1)*42+36 (i-1)*42+37 (i-1)*42+24;
        (i-1)*42+37 (i-1)*42+38 (i-1)*42+26;
        (i-1)*42+38 (i-1)*42+39 (i-1)*42+26;
        (i-1)*42+39 (i-1)*42+40 (i-1)*42+28;
        (i-1)*42+40 (i-1)*42+41 (i-1)*42+28;
        (i-1)*42+41 (i-1)*42+42 (i-1)*42+30;
        (i-1)*42+42 (i-1)*42+31 (i-1)*42+30;

        (i-1)*42+43 (i-1)*42+31 (i-1)*42+32;
        (i-1)*42+44 (i-1)*42+32 (i-1)*42+33;
        (i-1)*42+44 (i-1)*42+33 (i-1)*42+34;
        (i-1)*42+45 (i-1)*42+34 (i-1)*42+35;
        (i-1)*42+45 (i-1)*42+35 (i-1)*42+36;
        (i-1)*42+46 (i-1)*42+36 (i-1)*42+37;
        (i-1)*42+46 (i-1)*42+37 (i-1)*42+38;
        (i-1)*42+47 (i-1)*42+38 (i-1)*42+39;
        (i-1)*42+47 (i-1)*42+39 (i-1)*42+40;
        (i-1)*42+48 (i-1)*42+40 (i-1)*42+41;
        (i-1)*42+48 (i-1)*42+41 (i-1)*42+42;
        (i-1)*42+43 (i-1)*42+42 (i-1)*42+31;

        (i-1)*42+43 (i-1)*42+44 (i-1)*42+32;
        (i-1)*42+44 (i-1)*42+45 (i-1)*42+34;
        (i-1)*42+45 (i-1)*42+46 (i-1)*42+36;
        (i-1)*42+46 (i-1)*42+47 (i-1)*42+38;
        (i-1)*42+47 (i-1)*42+48 (i-1)*42+40;
        (i-1)*42+48 (i-1)*42+43 (i-1)*42+42;

        ];

    tri_direction=[tri_direction;
        ones(42,1)];
end

triNum=size(tri_ijk,1);

% Organize the trinagle sequence for the right direction
for i=1:triNum
    x1=node.coordinates_mat(tri_ijk(i,1),:);
    x2=node.coordinates_mat(tri_ijk(i,2),:);
    x3=node.coordinates_mat(tri_ijk(i,3),:);
    v1=x2-x1;
    v2=x3-x1;
    normal=cross(v1,v2);
    if sign(normal(3))==-tri_direction(i)
    else
        temp=tri_ijk(i,2);
        tri_ijk(i,2)=tri_ijk(i,3);
        tri_ijk(i,3)=temp;
    end
end

%% Define Triangle
cst.node_ijk_mat=tri_ijk;
cstNum=size(cst.node_ijk_mat,1);
cst.E_vec=E*ones(cstNum,1);
cst.t_vec=t*ones(cstNum,1);
cst.v_vec=v*ones(cstNum,1);

plots.Plot_Shape_CSTNumber;

%% Define Rotational Spring
for i=1:sectionNum
    rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
        (i-1)*42+18 (i-1)*42+7 (i-1)*42+1 (i-1)*42+8;
        (i-1)*42+7 (i-1)*42+8 (i-1)*42+1 (i-1)*42+2;
        (i-1)*42+1 (i-1)*42+8 (i-1)*42+2 (i-1)*42+9;
        (i-1)*42+8 (i-1)*42+9 (i-1)*42+2 (i-1)*42+10;
        (i-1)*42+9 (i-1)*42+10 (i-1)*42+2 (i-1)*42+3;
        (i-1)*42+2 (i-1)*42+10 (i-1)*42+3 (i-1)*42+11;
        (i-1)*42+10 (i-1)*42+11 (i-1)*42+3 (i-1)*42+12;
        (i-1)*42+11 (i-1)*42+12 (i-1)*42+3 (i-1)*42+4;
        (i-1)*42+3 (i-1)*42+12 (i-1)*42+4 (i-1)*42+13;
        (i-1)*42+12 (i-1)*42+13 (i-1)*42+4 (i-1)*42+14;
        (i-1)*42+13 (i-1)*42+14 (i-1)*42+4 (i-1)*42+5;
        (i-1)*42+4 (i-1)*42+14 (i-1)*42+5 (i-1)*42+15;
        (i-1)*42+14 (i-1)*42+15 (i-1)*42+5 (i-1)*42+16;
        (i-1)*42+15 (i-1)*42+16 (i-1)*42+5 (i-1)*42+6;
        (i-1)*42+5 (i-1)*42+16 (i-1)*42+6 (i-1)*42+17;
        (i-1)*42+16 (i-1)*42+17 (i-1)*42+6 (i-1)*42+18;
        (i-1)*42+17 (i-1)*42+18 (i-1)*42+6 (i-1)*42+1;
        (i-1)*42+6 (i-1)*42+18 (i-1)*42+1 (i-1)*42+7; 
        
        (i-1)*42+1 (i-1)*42+7 (i-1)*42+8 (i-1)*42+20;
        (i-1)*42+2 (i-1)*42+8 (i-1)*42+9 (i-1)*42+20;
        (i-1)*42+2 (i-1)*42+9 (i-1)*42+10 (i-1)*42+22;
        (i-1)*42+3 (i-1)*42+10 (i-1)*42+11 (i-1)*42+22;
        (i-1)*42+3 (i-1)*42+11 (i-1)*42+12 (i-1)*42+24;
        (i-1)*42+4 (i-1)*42+12 (i-1)*42+13 (i-1)*42+24;
        (i-1)*42+4 (i-1)*42+13 (i-1)*42+14 (i-1)*42+26;
        (i-1)*42+5 (i-1)*42+14 (i-1)*42+15 (i-1)*42+26;
        (i-1)*42+5 (i-1)*42+15 (i-1)*42+16 (i-1)*42+28;
        (i-1)*42+6 (i-1)*42+16 (i-1)*42+17 (i-1)*42+28;
        (i-1)*42+6 (i-1)*42+17 (i-1)*42+18 (i-1)*42+30;
        (i-1)*42+1 (i-1)*42+18 (i-1)*42+7 (i-1)*42+30;
        
        (i-1)*42+30 (i-1)*42+7 (i-1)*42+19 (i-1)*42+20;
        (i-1)*42+7 (i-1)*42+8 (i-1)*42+20 (i-1)*42+9;
        (i-1)*42+20 (i-1)*42+9 (i-1)*42+21 (i-1)*42+22;
        (i-1)*42+9 (i-1)*42+10 (i-1)*42+22 (i-1)*42+11;
        (i-1)*42+22 (i-1)*42+11 (i-1)*42+23 (i-1)*42+24;
        (i-1)*42+11 (i-1)*42+12 (i-1)*42+24 (i-1)*42+13;
        (i-1)*42+24 (i-1)*42+13 (i-1)*42+25 (i-1)*42+26;
        (i-1)*42+13 (i-1)*42+14 (i-1)*42+26 (i-1)*42+15;
        (i-1)*42+26 (i-1)*42+15 (i-1)*42+27 (i-1)*42+28;
        (i-1)*42+15 (i-1)*42+16 (i-1)*42+28 (i-1)*42+17;
        (i-1)*42+28 (i-1)*42+17 (i-1)*42+29 (i-1)*42+30;
        (i-1)*42+17 (i-1)*42+18 (i-1)*42+30 (i-1)*42+7;% Down

        (i-1)*42+42 (i-1)*42+31 (i-1)*42+43 (i-1)*42+32;
        (i-1)*42+31 (i-1)*42+32 (i-1)*42+43 (i-1)*42+44;
        (i-1)*42+43 (i-1)*42+32 (i-1)*42+44 (i-1)*42+33;
        (i-1)*42+32 (i-1)*42+33 (i-1)*42+44 (i-1)*42+34;
        (i-1)*42+33 (i-1)*42+34 (i-1)*42+44 (i-1)*42+45;
        (i-1)*42+44 (i-1)*42+34 (i-1)*42+45 (i-1)*42+35;
        (i-1)*42+34 (i-1)*42+35 (i-1)*42+45 (i-1)*42+36;
        (i-1)*42+35 (i-1)*42+36 (i-1)*42+45 (i-1)*42+46;
        (i-1)*42+45 (i-1)*42+36 (i-1)*42+46 (i-1)*42+37;
        (i-1)*42+36 (i-1)*42+37 (i-1)*42+46 (i-1)*42+38;
        (i-1)*42+37 (i-1)*42+38 (i-1)*42+46 (i-1)*42+47;
        (i-1)*42+46 (i-1)*42+38 (i-1)*42+47 (i-1)*42+39;
        (i-1)*42+38 (i-1)*42+39 (i-1)*42+47 (i-1)*42+40;
        (i-1)*42+39 (i-1)*42+40 (i-1)*42+47 (i-1)*42+48;
        (i-1)*42+47 (i-1)*42+40 (i-1)*42+48 (i-1)*42+41;
        (i-1)*42+40 (i-1)*42+41 (i-1)*42+48 (i-1)*42+42;
        (i-1)*42+41 (i-1)*42+42 (i-1)*42+48 (i-1)*42+43;
        (i-1)*42+48 (i-1)*42+42 (i-1)*42+43 (i-1)*42+31; 
        
        (i-1)*42+43 (i-1)*42+31 (i-1)*42+32 (i-1)*42+20;
        (i-1)*42+44 (i-1)*42+32 (i-1)*42+33 (i-1)*42+20;
        (i-1)*42+44 (i-1)*42+33 (i-1)*42+34 (i-1)*42+22;
        (i-1)*42+45 (i-1)*42+34 (i-1)*42+35 (i-1)*42+22;
        (i-1)*42+45 (i-1)*42+35 (i-1)*42+36 (i-1)*42+24;
        (i-1)*42+46 (i-1)*42+36 (i-1)*42+37 (i-1)*42+24;
        (i-1)*42+46 (i-1)*42+37 (i-1)*42+38 (i-1)*42+26;
        (i-1)*42+47 (i-1)*42+38 (i-1)*42+39 (i-1)*42+26;
        (i-1)*42+47 (i-1)*42+39 (i-1)*42+40 (i-1)*42+28;
        (i-1)*42+48 (i-1)*42+40 (i-1)*42+41 (i-1)*42+28;
        (i-1)*42+48 (i-1)*42+41 (i-1)*42+42 (i-1)*42+30;
        (i-1)*42+43 (i-1)*42+42 (i-1)*42+31 (i-1)*42+30;
        
        (i-1)*42+30 (i-1)*42+31 (i-1)*42+19 (i-1)*42+20;
        (i-1)*42+31 (i-1)*42+32 (i-1)*42+20 (i-1)*42+33;
        (i-1)*42+20 (i-1)*42+33 (i-1)*42+21 (i-1)*42+22;
        (i-1)*42+33 (i-1)*42+34 (i-1)*42+22 (i-1)*42+35;
        (i-1)*42+22 (i-1)*42+35 (i-1)*42+23 (i-1)*42+24;
        (i-1)*42+35 (i-1)*42+36 (i-1)*42+24 (i-1)*42+37;
        (i-1)*42+24 (i-1)*42+37 (i-1)*42+25 (i-1)*42+26;
        (i-1)*42+37 (i-1)*42+38 (i-1)*42+26 (i-1)*42+39;
        (i-1)*42+26 (i-1)*42+39 (i-1)*42+27 (i-1)*42+28;
        (i-1)*42+39 (i-1)*42+40 (i-1)*42+28 (i-1)*42+41;
        (i-1)*42+28 (i-1)*42+41 (i-1)*42+29 (i-1)*42+30;
        (i-1)*42+41 (i-1)*42+42 (i-1)*42+30 (i-1)*42+31; % Up
        ];

        if i<sectionNum
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
            (i-1)*42+32 (i-1)*42+43 (i-1)*42+44 (i-1)*42+50;
            (i-1)*42+34 (i-1)*42+44 (i-1)*42+45 (i-1)*42+52;
            (i-1)*42+36 (i-1)*42+45 (i-1)*42+46 (i-1)*42+54;
            (i-1)*42+38 (i-1)*42+46 (i-1)*42+47 (i-1)*42+56;
            (i-1)*42+40 (i-1)*42+47 (i-1)*42+48 (i-1)*42+58;
            (i-1)*42+42 (i-1)*42+48 (i-1)*42+43 (i-1)*42+60; % Connection
            ];
        end
end

% Find the bending stiffness
sprNum=size(rotSpr.node_ijkl_mat,1);
rotSpr.rot_spr_K_vec=kspr*ones(sprNum,1);
rotSpr.theta1=0.02*pi;
rotSpr.theta2=1.98*pi;
plots.Plot_Shape_SprNumber;

%% create the robot object
robot = loadrobot('universalUR5e','DataFormat','row');

% A default home position of the URe5 robot
% The 6*1 vector are angles of each joint
q_home1 = [0 -90 90 0 0 0]' * pi/180;


% This command creates the inverse kinematics sovlver for 
% The URe5 robots
ik = inverseKinematics('RigidBodyTree',robot);
ik.SolverParameters.AllowRandomRestart = false;
ikWeights1 = [1 1 1 1 1 1];

%% Initialize the assembly 
assembly.Initialize_Assembly;


%% Inflation and gravity application

targetP=30000;
nodalG=0.02;

nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    1 1 1 1;    
    2 1 1 1;
    3 1 1 1;
    4 1 1 1;
    5 1 1 1;
    6 1 1 1;
];

h=figure;
filename='UR_Pneumatic.gif';
pauseTime=0.01;
hold on

step=140;

Uhis=[];
for k=1:step
    
    xcurrent=node.coordinates_mat+node.current_U_mat;

    nodeNum=size(node.coordinates_mat,1);
    fmat=SolvePressureForce(xcurrent,k/step*targetP,tri_ijk);
    fmat(:,3)=fmat(:,3)-nodalG*k/step;
    fpressure=[(1:nodeNum)',fmat];

    nr.load=fpressure;
    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=10^-6;

    Uhis(k,:,:)=squeeze(nr.Solve());

    orient1 = [0 pi 0];
    position1 =  mean(xcurrent(1:6,:));
    
    tgtPose1 = eul2tform(orient1,'ZYX');
    tgtPose1(1:3,4)=position1';

    config1 = ik('tool0',double(tgtPose1),ikWeights1',q_home1');

    if mod(k,10)==0 
        clf
        
        set(gcf, 'color', 'white');
        
        hold on

        daspect([1 1 1]);
        xlim([0 0.8])
        ylim([-0.3 0.3])
        zlim([0.4 1])
        %view(30, 15)
        show(robot,config1)
        light("Style","local","Position",[10 -10 10]);
        
        plots.viewAngle1=0;
        plots.viewAngle2=0;

        plots.displayRange=[0 0.6 -0.3 0.3 0.4 1]';

        plots.Plot_DeformedShape(squeeze(Uhis(k,:,:)));

        pause(pauseTime);
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File 
        if k == 10 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 
        end
    end

end

UextEnd=squeeze(Uhis(end,:,:));
% figure
% plots.Plot_DeformedShape(UextEnd);

%% Rotation

sd=Solver_NR_SuppDisp;
sd.assembly=assembly;

sd.supp=[
    1 1 1 1;    
    2 1 1 1;
    3 1 1 1;
    4 1 1 1;
    5 1 1 1;
    6 1 1 1;
];

h=figure;
pauseTime=0.01;
hold on

step=140;

UhisRot=[];


xcurrent=node.coordinates_mat+node.current_U_mat;
Urot_current=node.current_U_mat;
baseNode_soft=xcurrent(1:6,:);
baseNodeCenter_soft=mean(baseNode_soft,1);
baseNodeRelative_soft=baseNode_soft-baseNodeCenter_soft;


for k=1:step

    % set up my target motion of soft robots
    sd.suppTarget = zeros(6, 4);

    angle_1=-0.002*pi*k; % Rotation by X-axis
    RotMat=[cos(-angle_1) 0 sin(-angle_1);
            0 1 0;   
            -sin(-angle_1) 0 cos(-angle_1)]; % Rotation by X-axis

    SuppU=baseNodeRelative_soft*RotMat-baseNodeRelative_soft;
    %baseNodeRelative_soft=baseNodeRelative_soft*RotMat;

    sd.suppTarget=[(1:6)', SuppU];
    
    xcurrent=node.coordinates_mat+node.current_U_mat;

    nodeNum=size(node.coordinates_mat,1);

    fmat=SolvePressureForce(xcurrent,targetP,tri_ijk);
    fmat(:,3)=fmat(:,3)-nodalG;
    node.current_ext_force_mat=fmat;
    
    sd.increStep=1;
    sd.iterMax=20;
    sd.tol=10^-6;
    
    Urot_current=squeeze(sd.Solve());
    % Uhis(k,:,:)=Urot_current+squeeze(sd.Solve());

    orient1 = [0 pi+angle_1 0];
    position1 =  mean(xcurrent(1:6,:));
    
    tgtPose1 = eul2tform(orient1,'ZYX');
    tgtPose1(1:3,4)=position1';

    config1 = ik('tool0',double(tgtPose1),ikWeights1',q_home1');

    if mod(k,10)==0 
        clf
        
        set(gcf, 'color', 'white');
        
        hold on

        daspect([1 1 1]);
        xlim([0 0.8])
        ylim([-0.3 0.3])
        zlim([0.4 1])
        %view(30, 15)
        show(robot,config1)
        light("Style","local","Position",[10 -10 10]);
        
        plots.viewAngle1=0;
        plots.viewAngle2=0;

        plots.displayRange=[0 0.6 -0.3 0.3 0.4 1]';

        % plots.Plot_DeformedShape(squeeze(Uhis(k,:,:)));
        plots.Plot_DeformedShape(Urot_current);

        pause(pauseTime);
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File 

        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 

    end

end

%% Bending Soft Structure

nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    1 1 1 1;
    2 1 1 1;
    3 1 1 1;
    4 1 1 1;
    5 1 1 1;
    6 1 1 1;
    19 1 1 1;
    20 1 1 1;
    21 1 1 1;
    22 1 1 1;
    23 1 1 1;
    24 1 1 1;
    25 1 1 1;
    26 1 1 1;
    27 1 1 1;
    28 1 1 1;
    29 1 1 1;
    30 1 1 1;];

step=100;
targetF_1=20; % Bending force on node 19 side
targetF_2=20;  % Bending force on node 21 side
targetF_3=0;  % Bending force on node 23 side
targetF_4=0;  % Bending force on node 25 side
targetF_5=0;  % Bending force on node 27 side
targetF_6=0;  % Bending force on node 29 side


frictionRate=0.8;
Uhis=[];

node.current_ext_force_mat=zeros(size(node.current_U_mat)); % Delete external force from last step
for k=1:step

    xcurrent=node.coordinates_mat+node.current_U_mat;
    nodeNum=size(node.coordinates_mat,1);
    fmat=SolvePressureForce(xcurrent,targetP,tri_ijk);
    fmat(:,3)=fmat(:,3)-nodalG;
    fpressure=[(1:nodeNum)',fmat];

    for i=1:sectionNum-1

        % Bending force on node 19 side
        n1=xcurrent((i)*42+19,:);
        n2=xcurrent((i-1)*42+19,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*42+19,2:4)=fpressure((i)*42+19,2:4)+v1*targetF_1*k/step*(frictionRate^(i));
        fpressure((i-1)*42+19,2:4)=fpressure((i-1)*42+19,2:4)+v2*targetF_1*k/step*(frictionRate^(i));

        % Bending force on node 21 side
        n1=xcurrent((i)*42+21,:);
        n2=xcurrent((i-1)*42+21,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*42+21,2:4)=fpressure((i)*42+21,2:4)+v1*targetF_2*k/step*(frictionRate^(i));
        fpressure((i-1)*42+21,2:4)=fpressure((i-1)*42+21,2:4)+v2*targetF_2*k/step*(frictionRate^(i));

        % Bending force on node 23 side
        n1=xcurrent((i)*42+23,:);
        n2=xcurrent((i-1)*42+23,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*42+23,2:4)=fpressure((i)*42+23,2:4)+v1*targetF_3*k/step*(frictionRate^(i));
        fpressure((i-1)*42+23,2:4)=fpressure((i-1)*42+23,2:4)+v2*targetF_3*k/step*(frictionRate^(i)); 

        % Bending force on node 25 side
        n1=xcurrent((i)*42+25,:);
        n2=xcurrent((i-1)*42+25,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*42+25,2:4)=fpressure((i)*42+25,2:4)+v1*targetF_4*k/step*(frictionRate^(i));
        fpressure((i-1)*42+25,2:4)=fpressure((i-1)*42+25,2:4)+v2*targetF_4*k/step*(frictionRate^(i));

        % Bending force on node 27 side
        n1=xcurrent((i)*42+27,:);
        n2=xcurrent((i-1)*42+27,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*42+27,2:4)=fpressure((i)*42+27,2:4)+v1*targetF_5*k/step*(frictionRate^(i));
        fpressure((i-1)*42+27,2:4)=fpressure((i-1)*42+27,2:4)+v2*targetF_5*k/step*(frictionRate^(i));

        % Bending force on node 29 side
        n1=xcurrent((i)*42+29,:);
        n2=xcurrent((i-1)*42+29,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*42+29,2:4)=fpressure((i)*42+29,2:4)+v1*targetF_6*k/step*(frictionRate^(i));
        fpressure((i-1)*42+29,2:4)=fpressure((i-1)*42+29,2:4)+v2*targetF_6*k/step*(frictionRate^(i));

     end

    nr.load=fpressure;
    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=10^-6;

    Ubend_current=squeeze(nr.Solve());

    orient1 = [0 pi+angle_1 0];
    position1 =  mean(xcurrent(1:6,:));
    
    tgtPose1 = eul2tform(orient1,'ZYX');
    tgtPose1(1:3,4)=position1';

    config1 = ik('tool0',double(tgtPose1),ikWeights1',q_home1');

    if mod(k,10)==0 
        clf
        
        set(gcf, 'color', 'white');
        
        hold on

        daspect([1 1 1]);
        xlim([0 0.8])
        ylim([-0.3 0.3])
        zlim([0.4 1])
        %view(30, 15)
        show(robot,config1)
        light("Style","local","Position",[10 -10 10]);
        
        plots.viewAngle1=0;
        plots.viewAngle2=0;

        plots.displayRange=[0 0.6 -0.3 0.3 0.4 1]';

        % plots.Plot_DeformedShape(squeeze(Uhis(k,:,:)));
        plots.Plot_DeformedShape(Ubend_current);

        pause(pauseTime);
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File 

        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 

    end

end



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