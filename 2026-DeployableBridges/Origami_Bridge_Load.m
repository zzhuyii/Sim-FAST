clear all
close all
clc
tic

%% Define Geometry, 
% Length of the bridge unit
L=2;

% Width of the bridge
W=4;

% Height of the bridge
H=2;

% Number of Sections
N=4;

% HSS 4X3X5/16 A500 Grade C Fy=50ksi
barA=0.0023; % 3.52 in^2
barE=2*10^11;

% We will use the weak axis
I=1.88*10^-6; 

% We assume a soft panel so that only the truss is taking global load
% Thus, panel Young's modulus is 200 MPa
panel_E=2*10^8;
panel_t=0.01;
panel_v=0.3;




%% Define assembly
node=Elements_Nodes;
assembly=Assembly_Origami;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;
bar=Vec_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_4N=rot_spr_4N;


%% Define Nodal Coordinates
for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        2*L*(i-1), 0, H;
        2*L*(i-1), 0, 0;
        2*L*(i-1), W, 0;
        2*L*(i-1), W, H;
        2*L*(i-1)+L, 0, H;
        2*L*(i-1)+L, 0, 0;
        2*L*(i-1)+L, L, 0;
        2*L*(i-1)+L, W, 0;
        2*L*(i-1)+L, W, H;];        
end

node.coordinates_mat=[node.coordinates_mat;
        2*L*N, 0, H;
        2*L*N, 0, 0;
        2*L*N, W, 0;
        2*L*N, W, H; 
        ];




%% Define Plotting Functions
plots=Plot_Kirigami_Truss;
plots.assembly=assembly;
plots.displayRange=[-1; 4*(N+1); -3; 7; -2; 3];

plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_Node_Number;


%% Define Triangle
for i=1:N
     cst.node_ijk_mat=[cst.node_ijk_mat;
        9*(i-1)+2    9*(i-1)+6    9*(i-1)+7;
        9*(i-1)+2    9*(i-1)+7    9*(i-1)+3;
        9*(i-1)+3    9*(i-1)+7    9*(i-1)+8;
        9*(i-1)+6    9*(i-1)+7    9*(i-1)+11;
        9*(i-1)+7    9*(i-1)+11    9*(i-1)+12;
        9*(i-1)+7    9*(i-1)+8    9*(i-1)+12;];
end

cstNum=size(cst.node_ijk_mat,1);
cst.t_vec=panel_t*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.v_vec=panel_v*ones(cstNum,1);

plots.Plot_Shape_CST_Number;



%% Define bar
for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        9*(i-1)+2   9*(i-1)+5;
        9*(i-1)+1   9*(i-1)+2;
        9*(i-1)+1   9*(i-1)+5;
        9*(i-1)+5   9*(i-1)+11;
        9*(i-1)+5   9*(i-1)+6;
        9*(i-1)+5   9*(i-1)+10;
        9*(i-1)+10   9*(i-1)+11;

        9*(i-1)+3   9*(i-1)+9;
        9*(i-1)+4   9*(i-1)+3;
        9*(i-1)+4   9*(i-1)+9;
        9*(i-1)+9   9*(i-1)+12;
        9*(i-1)+9   9*(i-1)+8;
        9*(i-1)+9   9*(i-1)+13;
        9*(i-1)+13   9*(i-1)+12;

        % add new bars
        9*(i-1)+2   9*(i-1)+7;
        9*(i-1)+3   9*(i-1)+7;
        9*(i-1)+6   9*(i-1)+7;
        9*(i-1)+8   9*(i-1)+7;
        9*(i-1)+11  9*(i-1)+7;
        9*(i-1)+12  9*(i-1)+7;
        9*(i-1)+2   9*(i-1)+3;
        9*(i-1)+3   9*(i-1)+8;
        9*(i-1)+8   9*(i-1)+12;
        9*(i-1)+2   9*(i-1)+6;
        9*(i-1)+6   9*(i-1)+11;
        ];
end

bar.node_ij_mat=[bar.node_ij_mat;
    9*(N-1)+11   9*(N-1)+12;
    ];

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);

plots.Plot_Shape_Bar_Number();

%% Define Rotational Spring
for i=1:N    
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        9*(i-1)+5  9*(i-1)+2  9*(i-1)+6  9*(i-1)+7;
        9*(i-1)+5  9*(i-1)+6  9*(i-1)+11  9*(i-1)+7;
        9*(i-1)+2  9*(i-1)+5  9*(i-1)+6  9*(i-1)+11;
        9*(i-1)+2  9*(i-1)+6  9*(i-1)+7  9*(i-1)+11;

        9*(i-1)+9  9*(i-1)+3  9*(i-1)+8  9*(i-1)+7;
        9*(i-1)+9  9*(i-1)+8  9*(i-1)+12  9*(i-1)+7;
        9*(i-1)+3  9*(i-1)+8  9*(i-1)+9  9*(i-1)+12;
        9*(i-1)+3  9*(i-1)+7  9*(i-1)+8  9*(i-1)+12;

        9*(i-1)+2  9*(i-1)+7  9*(i-1)+3  9*(i-1)+8;
        9*(i-1)+3  9*(i-1)+7  9*(i-1)+2  9*(i-1)+6;
        9*(i-1)+6  9*(i-1)+7  9*(i-1)+11  9*(i-1)+12;
        9*(i-1)+11  9*(i-1)+7  9*(i-1)+12  9*(i-1)+8;

        9*(i-1)+1  9*(i-1)+2  9*(i-1)+5  9*(i-1)+6;
        9*(i-1)+6  9*(i-1)+5  9*(i-1)+11  9*(i-1)+10;
        9*(i-1)+4  9*(i-1)+3  9*(i-1)+9  9*(i-1)+8;
        9*(i-1)+8  9*(i-1)+9  9*(i-1)+12  9*(i-1)+13;
        ];       
end

rotNumCenter=size(rot_spr_4N.node_ijkl_mat,1);

for i=1:N-1
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        9*(i-1)+5  9*(i-1)+10  9*(i-1)+11  9*(i-1)+14;
        9*(i-1)+9  9*(i-1)+12  9*(i-1)+13  9*(i-1)+18;
        9*(i-1)+7  9*(i-1)+11  9*(i-1)+12  9*(i-1)+16;
    ];        
end

rotNum=size(rot_spr_4N.node_ijkl_mat,1);
rot_spr_4N.rot_spr_K_vec=ones(rotNum,1)*10^8;

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;

assembly.Initialize_Assembly;


%% Calculate Self-weight of the Bridge

% Steel properties
rho_steel=7850;         % Density of steel in kg/m^3
g=9.81;                 % Gravitational acceleration (m/s^2)

% Cross-sectional area of bar elements
A_bar=barA;             % m^2 (consistent with bar.A_vec)

% Calculate the total length of all bar elements
L_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;

for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    p1=coords(n1,:);
    p2=coords(n2,:);
    len=norm(p1 - p2);
    L_total=L_total+len;
end

% Calculate the total weight of all bars (in Newtons)
W_bar=A_bar*L_total*rho_steel*g;



%% Set up solver  (Distributed load along full length on bridge bottom)
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[nr.supp;
     2      1 1 1;
     3      1 1 1;
     N*9+2  1 1 1;
     N*9+3  1 1 1;];


% force increment of each node per node
force=1000;   % N


for i=1:100

    nr.load=[];
    total_F=0;
    for k=1:N
        nr.load=[nr.load;
            6+(k-1)*9 0 0 -force*i;
            8+(k-1)*9 0 0 -force*i;
            11+(k-1)*9 0 0 -force*i;
            12+(k-1)*9 0 0 -force*i;
            ];
        total_F=force*4*i+total_F;
    end

    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=1e-5;
    
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
fprintf('Stiffness is: %.2f N/m\n', Kstiff);
fprintf('span/disp at failure is: %.2f \n', 16/Uaverage);
fprintf('capacity/weight: %.2f \n', total_F/W_bar);
fprintf('-----------------------------\n');

% Plot the bar stress
truss_stress=truss_strain.*(bar.E_vec);
plots.Plot_Shape_Bar_Stress(truss_stress);

% Plot failed bar stress
plots.Plot_Shape_Bar_Failure(passYN);



