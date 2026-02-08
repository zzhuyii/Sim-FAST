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
N=8;

% The cross section of this bridge is:
% HSS 4X3X5/16 A500 Grade C Fy=50ksi/200 GPa
barA=0.0023; 
barE=2*10^11;

% We will use the weak axis
I=1.88*10^-6; 

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
% Set up the plotting function for inspection
plots=Plot_Rolling_Bridge();
plots.assembly=assembly;

% We will plot for the Rolling Bridge
plots.displayRange=[-0.5;2*N+0.5;-0.5;2.5;-0.5;2.5]; 
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


%% Calculate Self-weight of the Bridge

% Steel properties
rho_steel=7850;         % Density of steel (kg/m^3)
g=9.81;                 % Gravitational acceleration (m/s^2)

% Bar cross-section area
A_bar=barA;             % m^2, should match bar.A_vec

% Total length of all normal bar elements
L_bar_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;

for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    p1=coords(n1,:);
    p2=coords(n2,:);
    len=norm(p1-p2);
    L_bar_total=L_bar_total+len;
end

% Total length of all actuator bar elements
L_actbar_total=0;
if isfield(actBar,'node_ij_mat')
    actBarNodeMat=actBar.node_ij_mat;
    for i=1:size(actBarNodeMat,1)
        n1=actBarNodeMat(i,1);
        n2=actBarNodeMat(i,2);
        p1=coords(n1,:);
        p2=coords(n2,:);
        len=norm(p1-p2);
        L_actbar_total=L_actbar_total+len;
    end
end

% Total bar length (all bars)
L_total=L_bar_total+L_actbar_total;

% Total weight of all bars (Newtons)
W_bar=A_bar*L_total*rho_steel*g;




%% Distributed load along full length on bridge bottom
nr=Solver_NR_Loading;
nr.assembly=assembly;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

% Set up the support of this bridge
nr.assembly=assembly;
nr.supp=[nodeNumVec,zeros(nodeNum,1),zeros(nodeNum,1),zeros(nodeNum,1)];
nr.supp(1,2:4)=ones(1,3);
nr.supp(4,2:4)=ones(1,3);
nr.supp(45,2:4)=ones(1,3);
nr.supp(47,2:4)=ones(1,3);

% force increment of each node per node
force=1000;   % N


for i=1:100

    % Nonlinear solver settings
    nr.increStep=1;
    nr.iterMax=50;
    nr.tol=1e-5;    

    nr.load=[];
    total_F=0;
    for k=1:N-1
        nr.load=[nr.load;
            6*(k-1)+3 0 0 -force*i;
            6*(k-1)+5 0 0 -force*i;];
        total_F=force*2*i+total_F;
    end
    
    % Solve
    Uhis=nr.Solve;
    

    %% Evaluate if member is failing
    % Deformation
    U_end=squeeze(Uhis(end,:,:));   % [nodeNum x 3]
    
    % strain in each bar
    truss_strain=bar.Solve_Strain(node, U_end); 

    % axial force in each bar (N)
    internal_force=truss_strain.*(bar.E_vec).*(bar.A_vec); 

    barNum=numel(internal_force);
    
    % effective length KL
    L0_vec=bar.L0_vec(:);
    K=1.0;                
    Lc=K.*L0_vec;
    
    % find r value: r = sqrt(I/A)
    r=sqrt(I/barA)*ones(barNum);
    % prevent division by zero
    r=max(r,1e-9); 
    
    % yield stress
    Fy=345*10^6;  % Q345 steel or Grade50 (345MPa,50ksi)
    
    % pass = 1 means member is not failing
    passYN=false(barNum,1);

    % Critical Stress Ratio
    StressRatio=NaN(barNum,1);

    % Failure mode
    modeStr=cell(barNum,1);

    % critical failure load
    Pn=NaN(barNum,1);

    % slenderness ratio: KL/r
    slender=NaN(barNum,1);   

    % Euler stress (Pa)
    Fe=NaN(barNum,1);   

    % Critical stress requirement from AISC (Pa)
    Fcr=NaN(barNum,1);   

    for k=1:barNum
        Ni=internal_force(k);
        Ai=bar.A_vec(k);
        Ei=bar.E_vec(k);
        Lci=Lc(k);
        ri=r(k);

        % Use the AISC to check the member failure
        [passi,modeStri,Pni,stressRatioi]=Check_Truss_AISC(Ni,Ai,Ei,Lci,ri,Fy);
    
        modeStr{k}=modeStri;
        Pn(k)=Pni;
        StressRatio(k)=stressRatioi;
        passYN(k)=passi;
    end

    if max(StressRatio)<1
        fprintf('All Truss Safe \n');
    else
        fprintf('Failure Detected \n');
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
fprintf('Failure load is: %.2f N\n', total_F);
fprintf('Mid-span deflection at failure is: %.3f m\n', Uaverage);
fprintf('Stiffness is: %.2f m\n', Kstiff);
fprintf('span/disp at failure is: %.2f \n', 16/Uaverage);
fprintf('capacity/weight: %.2f \n', total_F/W_bar);
fprintf('-----------------------------\n');

% Plot the bar stress
truss_stress=truss_strain.*(bar.E_vec);
plots.Plot_Shape_Bar_Stress(truss_stress);

% Plot failed bar stress
plots.Plot_Shape_Bar_Failure(passYN);

