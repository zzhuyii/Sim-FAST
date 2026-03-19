clear all;
clc;
close all;

%% Initialize the scissor, 
% Number of Sections
N=8;

% Height of the bridge
H=2; % (m)

% Length of the section
L=2; % (m)

% Primary Member
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
bt = 4/0.31;   % b/t ratio of wider face  (flange), from AISC tables
ht = 6/0.31;   % h/t ratio of narrower face (web),  from AISC tables

% Limiting slenderness for uniformly compressed plate elements in HSS
% AASHTO LRFD Art. 6.9.4.2.1, Table 6.9.4.2.1-1
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

% The three-node rotational spring stiffness
barL=sqrt(H^2+L^2);
kspr=barE*Iy/barL;



%% Define assembly
assembly=Assembly_Scissor_Bridge(); 
node=Elements_Nodes();
cst=Vec_Elements_CST;
rotSpr3N=CD_Elements_RotSprings_3N; 
rotSpr4N=Vec_Elements_RotSprings_4N;
bar=Vec_Elements_Bars;
actBar=Std_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_3N=rotSpr3N;
assembly.rot_spr_4N=rotSpr4N;
assembly.actBar=actBar;



%% Define the nodal coordinates
for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        L*(i-1), 0, 0;
        L*(i-1), L, 0;
        L*(i-1), 0, L;
        L*(i-1), L, L;

        (i-1)*L+L/2, 0, L/2;
        (i-1)*L+L/2, L, L/2;
        (i-1)*L+L/2, 0, 0;
        (i-1)*L+L/2, L, 0; 
           
        (i-1)*L+L/2, 0, L;
        (i-1)*L+L/2, L, L;
        ];
end

node.coordinates_mat=[node.coordinates_mat;
        L*N, 0, 0;
        L*N, L, 0;
        L*N, 0, L;
        L*N, L, L;
        ];


%% Define Plotting Functions
plots=Plot_Scissor_Bridge; 
plots.assembly=assembly;
plots.displayRange=[-1; 2*N+1; -1; 3; -1; 3]; 

plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_Node_Number;



%% Define Triangle
for i=1:N
    cst.node_ijk_mat=[cst.node_ijk_mat;
        10*(i-1)+1  10*(i-1)+2  10*(i-1)+7;
        10*(i-1)+2  10*(i-1)+7 10*(i-1)+8;
        10*(i-1)+7 10*(i-1)+8 10*(i-1)+11;
        10*(i-1)+8 10*(i-1)+12 10*(i-1)+11;
    ];
end

cstNum=size(cst.node_ijk_mat,1);
cst.t_vec=panel_t*ones(cstNum,1); 
cst.E_vec=panel_E*ones(cstNum,1); 
cst.v_vec=panel_v*ones(cstNum,1); 

plots.Plot_Shape_CST_Number;



%% Define bar
for i=1:N
    bar.node_ij_mat=[
        bar.node_ij_mat;
        10*(i-1)+1   10*(i-1)+7;
        10*(i-1)+7   10*(i-1)+11;
        10*(i-1)+2   10*(i-1)+8;
        10*(i-1)+8   10*(i-1)+12;

        10*(i-1)+3   10*(i-1)+9;
        10*(i-1)+9   10*(i-1)+13;
        10*(i-1)+4   10*(i-1)+10;
        10*(i-1)+10   10*(i-1)+14;

        10*(i-1)+3   10*(i-1)+4;
        10*(i-1)+9   10*(i-1)+10;

        10*(i-1)+1   10*(i-1)+2;
        10*(i-1)+7   10*(i-1)+8;

        10*(i-1)+1  10*(i-1)+5;
        10*(i-1)+3  10*(i-1)+5;
        10*(i-1)+2  10*(i-1)+6;
        10*(i-1)+4  10*(i-1)+6;  
        10*(i-1)+5  10*(i-1)+13;
        10*(i-1)+5  10*(i-1)+11;
        10*(i-1)+6  10*(i-1)+12;
        10*(i-1)+6  10*(i-1)+14;  

        10*(i-1)+3  10*(i-1)+10;
        10*(i-1)+4  10*(i-1)+9;
        10*(i-1)+13  10*(i-1)+10;
        10*(i-1)+9  10*(i-1)+14;

        10*(i-1)+2  10*(i-1)+7;
        10*(i-1)+1  10*(i-1)+8;
        10*(i-1)+8  10*(i-1)+11;
        10*(i-1)+7  10*(i-1)+12;
        ];

    bar.A_vec=[bar.A_vec;
        ones(8,1)*barA;
        ones(4,1)*barA_brace;
        ones(8,1)*barA;
        ones(8,1)*barA_brace;
        ];

end

i=N+1;
bar.node_ij_mat=[
    bar.node_ij_mat;
    10*(i-1)+1  10*(i-1)+2;
    10*(i-1)+3  10*(i-1)+4;     
    ];

bar.A_vec=[bar.A_vec;
    ones(2,1)*barA_brace;
    ];

barNum=size(bar.node_ij_mat);
barNum=barNum(1);
bar.E_vec=barE*ones(barNum,1);

plots.Plot_Shape_Node_Number();
plots.Plot_Shape_Bar_Number();



%% Define 3 Node Rotational Spring
for i=1:N    
    rotSpr3N.node_ijk_mat = [rotSpr3N.node_ijk_mat;
         10*(i-1)+1  10*(i-1)+5  10*(i-1)+13;
         10*(i-1)+3  10*(i-1)+5  10*(i-1)+11;
         10*(i-1)+2  10*(i-1)+6  10*(i-1)+14;
         10*(i-1)+4  10*(i-1)+6  10*(i-1)+12;

         10*(i-1)+4  10*(i-1)+3  10*(i-1)+5;
         10*(i-1)+3  10*(i-1)+4  10*(i-1)+6;
         10*(i-1)+2  10*(i-1)+1  10*(i-1)+5;
         10*(i-1)+6  10*(i-1)+2  10*(i-1)+1;

         10*(i-1)+5  10*(i-1)+11  10*(i-1)+12;
         10*(i-1)+11  10*(i-1)+12  10*(i-1)+6;
         10*(i-1)+5  10*(i-1)+13  10*(i-1)+14;
         10*(i-1)+13  10*(i-1)+14  10*(i-1)+6;

         10*(i-1)+11  10*(i-1)+13  10*(i-1)+14;
         10*(i-1)+13  10*(i-1)+14  10*(i-1)+12;
         10*(i-1)+14  10*(i-1)+12  10*(i-1)+11;
         10*(i-1)+12  10*(i-1)+11  10*(i-1)+13;

         ];
end

rotNum=size(rotSpr3N.node_ijk_mat,1);
rotSpr3N.rot_spr_K_vec=kspr*ones(rotNum,1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_RotSpr_3N_Number;



%% Set up four node rotational spring
for i=1:N    
    rotSpr4N.node_ijkl_mat = [rotSpr4N.node_ijkl_mat;
         10*(i-1)+1   10*(i-1)+2   10*(i-1)+7  10*(i-1)+8;
         10*(i-1)+7  10*(i-1)+8  10*(i-1)+11  10*(i-1)+12;

         10*(i-1)+4   10*(i-1)+3   10*(i-1)+10  10*(i-1)+9;
         10*(i-1)+10  10*(i-1)+9  10*(i-1)+14  10*(i-1)+13;
         ];
end

rotNum4N=size(rotSpr4N.node_ijkl_mat,1);
% Four-node rotational spring just for slab stability
rotSpr4N.rot_spr_K_vec=100000*ones(rotNum4N,1); 

plots.Plot_Shape_RotSpr_4N_Number;




%% Set up four node rotational spring
actBar.node_ij_mat=[];
for i=1:N
    actBar.node_ij_mat=[
        actBar.node_ij_mat;
        10*(i-1)+1 10*(i-1)+3;
        10*(i-1)+2 10*(i-1)+4;
        ];
end

i=N+1;
actBar.node_ij_mat=[
    actBar.node_ij_mat;
    10*(i-1)+1 10*(i-1)+3;
    10*(i-1)+2 10*(i-1)+4;
    ];

actBarNum1=size(actBar.node_ij_mat,1);

for i=1:N
    actBar.node_ij_mat=[
        actBar.node_ij_mat;
        10*(i-1)+7 10*(i-1)+5;
        10*(i-1)+5 10*(i-1)+9;
        10*(i-1)+6 10*(i-1)+10;
        10*(i-1)+8 10*(i-1)+6;
        ];
end

actBarNum=size(actBar.node_ij_mat,1);
actBar.A_vec=barA*ones(actBarNum,1);
actBar.E_vec=barE*ones(actBarNum,1);

plots.Plot_Shape_ActBar_Number;



%% Initialize Assembly
assembly.Initialize_Assembly;


%% Calculate Self-weight of the Bridge

% Steel properties
rho_steel=7850;         % Density of steel in kg/m^3
g=9.81;                 % Gravitational acceleration in m/s^2

% Bar elements
A_bar=barA;             % Cross-sectional area of bars in m^2 (matches bar.A_vec)

% Calculate total length of all bars
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

% AASHTO LRFD 3.6.1.6 pedestrian live load
qPL = 3.6e3;   % Pa 
W_LL = qPL*16*2; % N

% AASHTO Strength Combination
% 1.25 DC + 1.75 PL
W_factored=1.25*(W_bar+W_deck)+1.75*W_LL;



%% Set up solver  
nr=Solver_NR_Loading;
nr.assembly=assembly;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

nr.supp = [1    1 1 1;
           2    1 1 1;
           10*N+1    0 1 1; 
           10*N+2    0 1 1; 
           ];

% force increment of each node per node
force=W_factored/14/5;   % N


for i=1:5

    % Nonlinear solver settings
    nr.increStep=1;
    nr.iterMax=50;
    nr.tol=1e-5;  

    nr.load=[];
    total_F=0;
    for k=1:N-1
        nr.load=[nr.load;
            10*(k)+1 0 0 -force*i;
            10*(k)+2 0 0 -force*i;];
        total_F=force*2*i+total_F;
    end
    
    % Solve
    Uhis=nr.Solve;
    
    %% Evaluate if member is failing
    % Deformation
    U_end=squeeze(Uhis(end,:,:));   % [nodeNum x 3]
    
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
        fprintf('Step %2d : All Truss Members Safe (AASHTO LRFD)\n', i);
    else
        fprintf('Step %2d : Member Failure Detected (AASHTO LRFD)\n', i);
        break
    end

    % Check bending strength of scissor
    rotSpr3N.Solve_Global_Theta(node,U_end);
    thetaReal=rotSpr3N.theta_current_vec;
    thetaStressFree=rotSpr3N.theta_stress_free_vec;
    rotSpr3N_K=rotSpr3N.rot_spr_K_vec;

    MomentVec=abs(thetaReal-thetaStressFree).*rotSpr3N_K;
    maxMoment=1.5*max(MomentVec);
    % Dead load factor of 1.5;
    momentCapacity=Fy*Iy/0.0762;

    % -------------------------------------------------------------------
    % Mn= Rf Rb Rpc Myce
    % C6.12.2.2.2c
    % Rf Rb Rpc may be taken as 1.0 and Myce may be the fundamental yield
    % moment of the gross cross section. We have a doubly symmetric HSS
    % steel cross section. These shall be applied. 
    % -------------------------------------------------------------------

    if momentCapacity>maxMoment
        fprintf('Step %2d : Max Moment %2d < Capacity  %2d \n', i, maxMoment,momentCapacity);
    else 
        fprintf('Step %2d : Max Moment %2d > Capacity  %2d \n', i, maxMoment,momentCapacity);
        break
    end

end


% Find Stiffness
Uaverage=-mean(squeeze(Uhis(end,[3*N-3,3*N-1],3)));
Kstiff=total_F/Uaverage;

% Output results
fprintf('-----------------------------\n');
fprintf('Total length of all bars: %.2f m\n', L_total);
fprintf('Total bar weight: %.2f N\n', W_bar);
fprintf('Total load is: %.2f N\n', total_F);
fprintf('Mid-span deflection at Strength limit state is: %.3f m\n', Uaverage);
fprintf('Stiffness is: %.2f N/m\n', Kstiff);
fprintf('span/disp at Strength limit state is: %.2f \n', 16/Uaverage);
fprintf('Maximum DCR: %.2f \n', max(DCR));
fprintf('-----------------------------\n');

% Plot the bar stress
truss_stress=truss_strain.*(bar.E_vec);
plots.Plot_Shape_Bar_Stress(truss_stress);

% Plot failed bar stress
plots.Plot_Shape_Bar_Failure(passYN);

% Plot deformed shape
plots.Plot_Deformed_Shape(U_end);


