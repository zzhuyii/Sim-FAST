clear all;
clc;
close all;
tic

%% Define the Rolling Bridge Geometry
% Height of the bridge
H=2; % meter

% Width of the bridge
W=2; % meter

% Length of the section
L=2; % meter

% Number of Sections
N=8;

% Primary member
% HSS 8X4X5/16 A500 Grade C Fy=50ksi
barA=0.00415;
barE=2*10^11;

% Second moment of inertia in both direction
Iy=21.2*10^-6; 
Ix=7.16*10^-6;  

% Brace Member
% HSS 4X2X5/16
barA_brace=0.0019;

% -----------------------------------------------------------------------
%  AASHTO LRFD Material & Section Properties
%  (defined once outside the loop — these never change)
% -----------------------------------------------------------------------
Fy  = 345e6;    % Yield strength (Pa), ASTM A500 Gr.C
Fu  = 427e6;    % Tensile strength (Pa), A500 Gr.C
E   = barE;     % Elastic modulus (Pa)

% Resistance factors — AASHTO LRFD Art. 6.5.4.2
phi_ty = 0.95;  % Tension yielding (gross section)
phi_tf = 0.80;  % Tension fracture (net section)
phi_c  = 0.95;  % Axial compression

% Net section parameters (welded connection)
% AASHTO LRFD Art. 6.8.2.1
Ag = barA;
An = barA*0.9;   % Assume bolt hole is 0.1 of gross area
U  = 1.0;    % Shear lag factor
Rp = 1.0;    % Hole reduction factor (drilled/reamed holes) 

% Radius of gyration (weak axis governs)
r_val = sqrt(Ix/barA);

% Effective length factor
K = 1.0;

%  Section: HSS 4x3x5/16, A500 Gr.C
bt = 7.31;   % b/t ratio of wider face  (flange), from AISC tables
ht = 10.7;   % h/t ratio of narrower face (web),  from AISC tables

% Limiting slenderness for uniformly compressed plate elements in HSS
% AASHTO LRFD Art. 6.9.4.2.1, Table 6.9.4.2.1-1
% lambda_r = 1.92 * sqrt(E / Fy);
lambda_r = 1.28 * sqrt(E / Fy);

local_buckle_bt_pass = (bt <= lambda_r);
local_buckle_ht_pass = (ht <= lambda_r);
local_buckle_pass    = local_buckle_bt_pass && local_buckle_ht_pass;

fprintf('--- Local Buckling Check (AASHTO LRFD Art. 6.9.4.2) ---\n');
if local_buckle_pass
    fprintf('  Section is non-slender (local buckling OK)\n');
else
    fprintf('  WARNING: Section FAILS local buckling slenderness limit\n');
end


% We assume a soft panel so that only the truss is taking global load
% Thus, panel Young's modulus is 200 MPa
panel_E=2*10^8;
panel_t=0.01;
panel_v=0.3;

% We assume that the actuator bar has similar material property
% and stiffness as normal truss
activeBarE=2*10^11;

% Load the deformation history
UhisNew=load('RollingUhis.mat');
UhisNew=UhisNew.Uhis;
DepRate=0.3; % 1 is fully deployed, 0 is compact
DepStep = max(1, int32((1 - DepRate) * 800));


%% Initialize Elements and Assembly
node=Elements_Nodes;
bar=Vec_Elements_Bars;
actBar=CD_Elements_Bars;
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


%% Define Plotting Functions
plots=Plot_Rolling_Bridge();
plots.assembly=assembly;

plots.displayRange=[-2;18;-1;3;-1;14];
plots.viewAngle1=20;
plots.viewAngle2=20;

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

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);

% Add top bars for stability
bar.node_ij_mat=[bar.node_ij_mat;
    2 6;
    7 8;
    ];

bar.A_vec=[bar.A_vec;
    ones(2,1)*barA_brace;
    ];

for i=1:N-1
    bar.node_ij_mat=[bar.node_ij_mat;
        10+(i-2)*6,  12+(i-2)*6;
        13+(i-2)*6,  14+(i-2)*6;
        16+(i-2)*6,  18+(i-2)*6;
        ];
    bar.A_vec=[bar.A_vec;
        ones(3,1)*barA_brace;
        ];
end

barNum=size(bar.node_ij_mat,1);
bar.E_vec=barE*ones(barNum,1);

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
node.coordinates_mat=node.coordinates_mat+squeeze(UhisNew(DepStep,:,:));
assembly.Initialize_Assembly();


%% Calculate Self-weight of the Bridge

% Steel properties
rho_steel=7850;         % Density of steel (kg/m^3)
g=9.81;                 % Gravitational acceleration (m/s^2)

% Total length of all normal bar elements
L_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;

% Total bar weight (Newtons)
W_bar=0; 

for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    p1=coords(n1,:);
    p2=coords(n2,:);
    len=norm(p1-p2);

    L_total=L_total+len;
    W_bar=W_bar+(len*bar.A_vec(i)*rho_steel*g);
end

% -------------------------------------------------------------------
% Assume the deck is made with a 3 cm thick wood plate
% Supported using a 20 cm by 10 cm beam with 50 cm spacing
% The total weight of deck is:
% -------------------------------------------------------------------
W_deck=2*(0.03+10/50*0.2)*16*1000*9.8;


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

% Set up the support of this bridge
nr.assembly=assembly;
nr.supp=[nodeNumVec,zeros(nodeNum,1),zeros(nodeNum,1),zeros(nodeNum,1)];
nr.supp(1,2:4)=ones(1,3);
nr.supp(4,2:4)=ones(1,3);
nr.supp(3,2:4)=ones(1,3);
nr.supp(5,2:4)=ones(1,3);




for i = 1:5

    % Nonlinear solver settings
    nr.increStep = 1;
    nr.iterMax   = 50;
    nr.tol       = 1e-5;

    % Apply self weight incrementally
    nodeNum = size(node.coordinates_mat, 1);
    force   = (W_bar+W_deck) / nodeNum / 5 * i;

    nr.load = [(1:nodeNum)'  zeros(nodeNum,1) ...
                zeros(nodeNum,1)  -force*ones(nodeNum,1)];

    % Solve
    Uhis = nr.Solve;


    %%  Evaluate if member is failing (AASHTO LRFD)
    % Deformation
    U_end = squeeze(Uhis(end,:,:));   % [nodeNum x 3]

    % Strain and internal force in each bar
    truss_strain   = bar.Solve_Strain(node, U_end);
    internal_force = truss_strain .* (bar.E_vec) .* (bar.A_vec);

    % effective length KL
    L0_vec  = bar.L0_vec(:);
    Lc      = K .* L0_vec;

    % Output arrays (reset each load step)
    barNum=numel(internal_force);
    passYN      = false(barNum, 1);  % pass = 1 means member is not failing
    DCR         = NaN(barNum, 1);  % Demand over capacity ratio (Critical Stress) Ratio
    phi_Rn      = NaN(barNum, 1);  % Failure mode
    modeStr     = cell(barNum, 1);  % Failure mode
    slender_chk = cell(barNum, 1);  % Slenderness ratio check

    % Member-by-member check
    for k = 1:barNum

        Pu_k=1.5*internal_force(k);
        % Dead load factor use 1.5
        Ag_k=bar.A_vec(k);
        An_k=An;
        E_k=bar.E_vec(k);
        KL_k=Lc(k);
        r_k=r_val;
        Fy_k=Fy;
        Fu_k=Fu;
        Rp_k=Rp;      

        % Use the AASTHO to check the member failure
        [pass_k,modeStr_k,Pn_k,phi,phiPn,DCR_k]=Check_Truss_LRFD(Pu_k,...
            Ag_k, An_k, E_k, KL_k, r_k, Fy_k, Fu_k, Rp_k);
    
        passYN(k)=pass_k;
        modeStr{k}=modeStr_k;
        Pn(k)=Pn_k;
        DCR(k)=DCR_k; 

    end  

    if all(passYN)
        fprintf('Step %2d/50 : All Truss Members Safe (AASHTO LRFD)\n', i);
    else
        fprintf('Step %2d/50 : Member Failure Detected (AASHTO LRFD)\n', i);
        break
    end

end  

% Find Stiffness
Uaverage=-mean(squeeze(Uhis(end,[45,47],3)));

% Output results
fprintf('-----------------------------\n');
fprintf('Total length of all bars: %.2f m\n', L_total);
fprintf('Total bar weight: %.2f N\n', W_bar);
fprintf('Maximum stress ratio: %.2f \n', max(DCR));
fprintf('Tip deflection: %.2f \n', Uaverage);
fprintf('-----------------------------\n');



% Plot bar stress
truss_stress = truss_strain .* (bar.E_vec);
plots.Plot_Shape_Bar_Stress(truss_stress);

% Plot failed bars
plots.Plot_Shape_Bar_Failure(passYN);

% Plot deformed shape
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)));


