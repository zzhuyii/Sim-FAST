clear all
close all
clc

%% Define Geometry
sectionNum=5;

% We make the hexgonal robot to have tapered cross section
Bottom_Hex_side_length=0.05;
side_shrink_rate=0.95;
Hex_side_length=[];
for i=1:sectionNum
Hex_side_length=[Hex_side_length
         Bottom_Hex_side_length*(side_shrink_rate^(i-1));];
end

% Other geometries
Hole_hex_side_length=0.01; % Edge length of pneumatic hole at the center 
layerH=0.002; % Layer height initially

% film properties
t=0.00011; % 0.11mm for the thickness of film
E=1.3*10^9; % Young's modulus of the film
v=0.2; % Poisson's Ratio of the film

% rotational spring's stiffness
kspr=10*E*t^3/12;

% target pressure
targetP=40000;

%% Define the Nodal Coordinates
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
plots.viewAngle1=70;
plots.viewAngle2=15;
plots.displayRange=[-0.15; 0.15; -0.15; 0.15; -0.02; 0.3];

% plot the nodal location for inspection
plots.Plot_Shape_NodeNumber;

%% Set up the Triangle Elements
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

    % The bottom film is facing down
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

    % The top film is facing up
    tri_direction=[tri_direction;
        ones(42,1)];
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

% Plot the rotational springs for inspection
plots.Plot_Shape_SprNumber;

%% Initialize the assembly
assembly.Initialize_Assembly;


%% Set up solver for inflation
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

step=50;
Uhis=[];

for k=1:step

    % To help with convergence, we further make the first few steps 
    % to have smaller step lengths 
    if k<2
        subStep=20;
        for q=1:subStep
            % Update the pressure forces for nonlinearity
            xcurrent=node.coordinates_mat+node.current_U_mat;    
            nodeNum=size(node.coordinates_mat,1);
            fmat=SolvePressureForce(xcurrent,(k-1+q/subStep)/step*targetP,tri_ijk);
            fpressure=[(1:nodeNum)',fmat];
        
            nr.load=fpressure;            
            nr.increStep=1;
            nr.iterMax=20;
            nr.tol=10^-8;
        
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


%% Set up solver for the bending experiments
nr.supp=[
    19 1 1 1;
    21 1 1 1;
    23 1 1 1;
    25 1 1 1;
    27 1 1 1;
    29 1 1 1;
];

step=100;

frictionRate=0.86;

targetF=0;
deltaF=0.1;

UhisBend=[];

for p=1:10

    targetBend = input('Please enter a value for target bending angle (degree, between 0 and 50): '); 

    for k=1:step

        xcurrent=node.coordinates_mat+node.current_U_mat;

        node1=xcurrent(213,:);
        node2=xcurrent(215,:);
   
        normal=node2-node1;
        normal=normal/norm(normal);

        zaxis=[0 0 1];       

        currentBendAngle=(acos(dot(zaxis,normal))-pi/2)/pi*180;

        if currentBendAngle<targetBend
            targetF=targetF+deltaF;
        else
            targetF=targetF-deltaF;
        end
    
        % Update the pressure forces for nonlinearity
        
        nodeNum=size(node.coordinates_mat,1);
        fmat=SolvePressureForce(xcurrent,targetP,tri_ijk);
        fpressure=[(1:nodeNum)',fmat];  
    
        % Apply the cable tension as forces
        for i=1:sectionNum-1
            n1=xcurrent((i)*42+27,:);
            n2=xcurrent((i-1)*42+27,:);
    
            v1=n2-n1;
            v2=n1-n2;
    
            v1=v1/norm(v1);
            v2=v2/norm(v2);
    
            fpressure((i)*42+27,2:4)=fpressure((i)*42+27,2:4)+v1*targetF;
            fpressure((i-1)*42+27,2:4)=fpressure((i-1)*42+27,2:4)+v2*targetF; % Bending force on node 27 side
    
            n1=xcurrent((i)*42+29,:);
            n2=xcurrent((i-1)*42+29,:);
    
            v1=n2-n1;
            v2=n1-n2;
    
            v1=v1/norm(v1);
            v2=v2/norm(v2);
    
            fpressure((i)*42+29,2:4)=fpressure((i)*42+29,2:4)+v1*targetF;
            fpressure((i-1)*42+29,2:4)=fpressure((i-1)*42+29,2:4)+v2*targetF; % Bending force on node 29 side    
    
         end
    
        nr.load=fpressure;    
        nr.increStep=1;
        nr.iterMax=20;
        nr.tol=10^-6;
    
        UhisBend(k,:,:)=squeeze(nr.Solve());        
    
    end    

    fprintf('To achieve target bending angle of %d degree \n', targetBend)
    fprintf('The target cable tension is: %s N \n',  targetF)
    figure
    plots.Plot_DeformedShape(squeeze(UhisBend(end,:,:)))

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


