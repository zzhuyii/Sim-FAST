clear all;
clc;
close all;

%% Define the Rolling Bridge Geometry
% Height of the bridge
H=2; % meter

% Width of the bridge
W=2; % meter

% Length of the section
L=2; % meter

% Number of Sections
N=2;

% The cross section of this bridge is:
% HSS 8X4X5/16 A500 Grade C Fy=50ksi
barA=0.00415;
barE=2*10^11;

% We assume a soft panel so that only the truss is taking global load
% Thus, panel Young's modulus is 200 MPa
panel_E=2*10^8;
panel_t=0.01;
panel_v=0.3;

% We assume that the actuator bar has similar material property
% and stiffness as normal truss
activeBarE=2*10^11; % 80*10^9;

% Initialize Elements and Assembly
node=Elements_Nodes;
bar=Vec_Elements_Bars;
actBar=Std_Elements_Bars;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;

assembly=Assembly_Rolling_Bridge;
assembly.node=node;
assembly.bar=bar;
assembly.actBar=actBar;
assembly.cst=cst;
assembly.rot_spr_4N=rot_spr_4N;



%% Define the nodal coordinates
node.coordinates_mat=[];
node.coordinates_mat=[
    node.coordinates_mat;
    0    0   0;
    L/2  0   H;
    L    0   0;
    0    W   0;
    L    W   0;
    L/2  W   H;
    L    0   H;
    L    W   H;
    ];

for i=2:N-1
    node.coordinates_mat=[node.coordinates_mat;
        L*i      0   0;
        L*i-L/2  0   H;
        L*i      W   0;
        L*i-L/2  W   H;
        L*i      0   H;
        L*i      W   H;
        ];
end

node.coordinates_mat=[node.coordinates_mat;
    L*N      0   0;
    L*N-L/2  0   H;
    L*N      W   0;
    L*N-L/2  W   H;
    ];

% Set up the plotting function for inspection
plots=Plot_Rolling_Bridge();
plots.assembly=assembly;

% We will plot for the Rolling Bridge
plots.displayRange=[-2;18;-1;3;-1;14]; 
plots.viewAngle1=20;
plots.viewAngle2=20;

% Plot the nodal coordinates for inspection
plots.Plot_Shape_Node_Number()


%% Define triangle elements
cst.node_ijk_mat=[cst.node_ijk_mat;
    1  3  4;
    3  4  5;
    ];

for i=2:N
    cst.node_ijk_mat=[cst.node_ijk_mat;
        3+(i-2)*6  5+(i-2)*6  9+(i-2)*6;
        5+(i-2)*6  9+(i-2)*6  11+(i-2)*6;
        ];
end    

cstNum=size(cst.node_ijk_mat,1);
cst.v_vec=panel_v*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.t_vec=panel_t*ones(cstNum,1);

% Plot triangle elements for inspection
plots.Plot_Shape_CST_Number();


%% Define how normal bars are connected
bar.node_ij_mat=[bar.node_ij_mat;
    1, 2;
    1, 3;
    2, 3;
    4, 5;
    4, 6;
    5, 6;
    2, 7;
    6, 8;
    1, 4;
    3, 5;
    3, 4;
    ];

for i=2:N-1
    bar.node_ij_mat=[bar.node_ij_mat;
        3+(i-2)*6,  9+(i-2)*6;
        3+(i-2)*6,  10+(i-2)*6;
        9+(i-2)*6,  10+(i-2)*6;
        5+(i-2)*6,  11+(i-2)*6;
        5+(i-2)*6,  12+(i-2)*6;
        11+(i-2)*6, 12+(i-2)*6;
        7+(i-2)*6,  10+(i-2)*6;
        8+(i-2)*6,  12+(i-2)*6;
        10+(i-2)*6, 13+(i-2)*6;
        12+(i-2)*6, 14+(i-2)*6;
        9+(i-2)*6,  11+(i-2)*6;
        5+(i-2)*6,  9+(i-2)*6;
        ];
end

bar.node_ij_mat=[bar.node_ij_mat;
    3+(N-2)*6,  9+(N-2)*6;
    3+(N-2)*6,  10+(N-2)*6;
    9+(N-2)*6,  10+(N-2)*6;
    5+(N-2)*6,  11+(N-2)*6;
    5+(N-2)*6,  12+(N-2)*6;
    11+(N-2)*6, 12+(N-2)*6;
    7+(N-2)*6,  10+(N-2)*6;
    8+(N-2)*6,  12+(N-2)*6;
    5+(N-2)*6,  9+(N-2)*6;
    9+(N-2)*6,  11+(N-2)*6;
    ];

% Add top bars for stability
bar.node_ij_mat=[bar.node_ij_mat;
    2 6;
    7 8;
    ];

for i=1:N-1
    bar.node_ij_mat=[bar.node_ij_mat;
        10+(i-2)*6,  12+(i-2)*6;
        13+(i-2)*6,  14+(i-2)*6;
        16+(i-2)*6,  18+(i-2)*6;
        ];
end

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);

% Plot bar elements for inspection
plots.Plot_Shape_Bar_Number();


%% Define how actuator bars are connected
for i=1:N-1
    actBar.node_ij_mat=[actBar.node_ij_mat;
        3+(i-1)*6, 7+(i-1)*6;
        5+(i-1)*6, 8+(i-1)*6;
        ];
end

actBarNum=size(actBar.node_ij_mat,1);
actBar.A_vec=barA*ones(actBarNum,1);
actBar.E_vec=activeBarE*ones(actBarNum,1);

% Plot actuator bar elements for inspection
plots.Plot_Shape_ActBar_Number();


%% Define four node rotational springs

% Major Rot Spr
rot_spr_4N.node_ijkl_mat=[
    4 1 3 2;
    3 4 5 6;];

for i=1:N-1
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        6*(i-1)+5, 6*(i-1)+3, 6*(i-1)+9, 6*(i-1)+10;
        6*(i-1)+9, 6*(i-1)+5, 6*(i-1)+11, 6*(i-1)+12;];
end

% Floor Rot Spr
rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    1 3 4 5];

for i=1:N-1
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        6*(i-1)+3, 6*(i-1)+9, 6*(i-1)+5, 6*(i-1)+11;];
end

% Secondary Rot Spring
rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    1 2 3 7;
    4 5 6 8;];

for i=1:N-2
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        6*(i-1)+7, 6*(i-1)+3, 6*(i-1)+10, 6*(i-1)+9;
        6*(i-1)+3, 6*(i-1)+9, 6*(i-1)+10, 6*(i-1)+13;
        6*(i-1)+8, 6*(i-1)+5, 6*(i-1)+12, 6*(i-1)+11;
        6*(i-1)+5, 6*(i-1)+11, 6*(i-1)+12, 6*(i-1)+14;];
end

i=N-1;
rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    6*(i-1)+7, 6*(i-1)+3, 6*(i-1)+10, 6*(i-1)+9;
    6*(i-1)+8, 6*(i-1)+5, 6*(i-1)+12, 6*(i-1)+11;
    ];

rotNum=size(rot_spr_4N.node_ijkl_mat,1);
rot_spr_4N.rot_spr_K_vec=(10^6)*ones(rotNum,1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;


%% Initialize the entire assembly 
assembly.Initialize_Assembly();


%% Set up the self actuation solver
ta=Solver_NR_TrussAction;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

% Set up the support of this bridge
% The left end of this bridge is fully restricted
ta.assembly=assembly;
ta.supp=[nodeNumVec,zeros(nodeNum,1),zeros(nodeNum,1),zeros(nodeNum,1)];
ta.supp(1,2:4)=ones(1,3);
ta.supp(2,2:4)=ones(1,3);
ta.supp(3,2:4)=ones(1,3);
ta.supp(4,2:4)=ones(1,3);

% Set up the total loading step
ta.increStep=800;
% Set up the maximum iteration
ta.iterMax=30;
% Set up the tolorence
ta.tol=10^-1;

% Extension of actuator bars
dL=1.1;
ta.targetL0=actBar.L0_vec;
ta.targetL0=ta.targetL0+dL;

% Solve for the deformation history
Uhis=ta.Solve();

% Plot the deformed shape
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)));

% Also plot the deformation history
plots.fileName="Rolling_Bridge_Deploy_Unit.gif";
plots.Plot_Deformed_His(Uhis(1:10:end,:,:));


% store the deformation history
% save('RollingUhis.mat','Uhis'); % Saves to a .mat file
% UhisNew=load('RollingUhis.mat');
% UhisNew=UhisNew.Uhis;