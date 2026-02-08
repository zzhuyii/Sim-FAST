clear all
close all
clc
tic

%% Define Geometry
% Length of the section
L=2;

% Gap for setting up bars
gap=0;

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



%% Define assembly
assembly=Assembly_Kirigami_Truss;
node=Elements_Nodes;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;
rot_spr_3N=CD_Elements_RotSprings_3N;
bar=Vec_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_4N=rot_spr_4N;
assembly.rot_spr_3N=rot_spr_3N;


%% Define Nodes

node.coordinates_mat=[node.coordinates_mat;
        0, 0, 0;
        0, L, 0;
        0, 0, L;
        0, L, L;];

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        (L)*(i-1)+L/2, 0, 0;
        (L)*(i-1)+L/2, 0, gap;
        (L)*(i-1)+L/2, L, 0;
        (L)*(i-1)+L/2, L, gap;

        (L)*(i-1)+L/2, 0, L;
        (L)*(i-1)+L/2, gap, L;
        (L)*(i-1)+L/2, L, L;
        (L)*(i-1)+L/2, L, L-gap;

        (L)*(i-1)+L/2, L/2, 0;
        (L)*(i-1)+L/2, L/2, L;
        (L)*(i-1)+L/2, 0, L/2;
        (L)*(i-1)+L/2, L, L/2;

        (L)*(i-1)+L, 0, 0;
        (L)*(i-1)+L, L, 0;
        (L)*(i-1)+L, 0, L;
        (L)*(i-1)+L, L, L;];
end



%% Define Plotting Functions
plots=Plot_Kirigami_Truss;
plots.assembly=assembly;
plots.displayRange=[-1; L*(N+1); -1; 3; -1; 3];

plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_Node_Number;


%% Define Triangle
for i=1:N
     cst.node_ijk_mat=[cst.node_ijk_mat;
        16*(i-1)+1    16*(i-1)+5    16*(i-1)+13;
        16*(i-1)+1    16*(i-1)+13    16*(i-1)+2;
        16*(i-1)+2    16*(i-1)+7    16*(i-1)+13;
        16*(i-1)+7    16*(i-1)+13    16*(i-1)+18;
        16*(i-1)+13    16*(i-1)+17    16*(i-1)+18;
        16*(i-1)+13    16*(i-1)+5    16*(i-1)+17;
        ];
end

cstNum=size(cst.node_ijk_mat,1);
cst.t_vec=panel_t*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.v_vec=panel_v*ones(cstNum,1);

plots.Plot_Shape_CST_Number;

%% Define bar
for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        16*(i-1)+1   16*(i-1)+2;
        16*(i-1)+1   16*(i-1)+3;
        16*(i-1)+2   16*(i-1)+4;
        16*(i-1)+3   16*(i-1)+4;

        16*(i-1)+1   16*(i-1)+5;
        16*(i-1)+5   16*(i-1)+17;
        16*(i-1)+2   16*(i-1)+7;
        16*(i-1)+7   16*(i-1)+18;
        
        16*(i-1)+3   16*(i-1)+9;
        16*(i-1)+9   16*(i-1)+19;
        16*(i-1)+4   16*(i-1)+11;
        16*(i-1)+11   16*(i-1)+20;

        16*(i-1)+4   16*(i-1)+12;
        16*(i-1)+12   16*(i-1)+20;
        16*(i-1)+3   16*(i-1)+10;
        16*(i-1)+10   16*(i-1)+19;

        16*(i-1)+1   16*(i-1)+6;
        16*(i-1)+6   16*(i-1)+17;
        16*(i-1)+2   16*(i-1)+8;
        16*(i-1)+8   16*(i-1)+18;

        16*(i-1)+1   16*(i-1)+15;
        16*(i-1)+3   16*(i-1)+15;
        16*(i-1)+15   16*(i-1)+17;
        16*(i-1)+15   16*(i-1)+19;

        16*(i-1)+2   16*(i-1)+16;
        16*(i-1)+16   16*(i-1)+20;
        16*(i-1)+4   16*(i-1)+16;
        16*(i-1)+16   16*(i-1)+18;

        16*(i-1)+3   16*(i-1)+14;
        16*(i-1)+4   16*(i-1)+14;
        16*(i-1)+14   16*(i-1)+20;
        16*(i-1)+14   16*(i-1)+19;

        16*(i-1)+1   16*(i-1)+13;
        16*(i-1)+2   16*(i-1)+13;
        16*(i-1)+13   16*(i-1)+18;
        16*(i-1)+13   16*(i-1)+17;

        16*(i-1)+10   16*(i-1)+14;
        16*(i-1)+11   16*(i-1)+14;
        16*(i-1)+12   16*(i-1)+16;
        16*(i-1)+8   16*(i-1)+16;

        16*(i-1)+7   16*(i-1)+13;
        16*(i-1)+5   16*(i-1)+13;
        16*(i-1)+6   16*(i-1)+15;
        16*(i-1)+15   16*(i-1)+9;
        ];
end

i=N+1;
bar.node_ij_mat=[bar.node_ij_mat;
    16*(i-1)+1   16*(i-1)+2;
    16*(i-1)+1   16*(i-1)+3;
    16*(i-1)+2   16*(i-1)+4;
    16*(i-1)+3   16*(i-1)+4;];

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);

plots.Plot_Shape_Bar_Number();
plots.Plot_Shape_Node_Number();


%% Define Rotational Spring
for i=1:N    
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        16*(i-1)+1  16*(i-1)+6  16*(i-1)+15  16*(i-1)+17;
        16*(i-1)+3  16*(i-1)+9  16*(i-1)+15  16*(i-1)+19;

        16*(i-1)+1  16*(i-1)+3  16*(i-1)+15  16*(i-1)+9;
        16*(i-1)+3  16*(i-1)+1  16*(i-1)+15  16*(i-1)+6;
        16*(i-1)+6  16*(i-1)+15  16*(i-1)+17  16*(i-1)+19;
        16*(i-1)+9  16*(i-1)+15  16*(i-1)+19  16*(i-1)+17;

        16*(i-1)+3  16*(i-1)+10  16*(i-1)+14  16*(i-1)+19;
        16*(i-1)+4  16*(i-1)+11  16*(i-1)+14  16*(i-1)+20;

        16*(i-1)+11  16*(i-1)+14  16*(i-1)+4  16*(i-1)+3;
        16*(i-1)+4  16*(i-1)+14  16*(i-1)+3  16*(i-1)+10;
        16*(i-1)+10  16*(i-1)+14  16*(i-1)+19  16*(i-1)+20;
        16*(i-1)+11  16*(i-1)+14  16*(i-1)+20  16*(i-1)+19;

        16*(i-1)+2  16*(i-1)+8  16*(i-1)+16  16*(i-1)+18;
        16*(i-1)+4  16*(i-1)+12  16*(i-1)+16  16*(i-1)+20;

        16*(i-1)+2  16*(i-1)+16  16*(i-1)+4  16*(i-1)+12;
        16*(i-1)+4  16*(i-1)+16  16*(i-1)+2  16*(i-1)+8;
        16*(i-1)+8  16*(i-1)+16  16*(i-1)+18  16*(i-1)+20;
        16*(i-1)+18  16*(i-1)+16  16*(i-1)+20  16*(i-1)+12;

        16*(i-1)+2  16*(i-1)+7  16*(i-1)+13  16*(i-1)+18;
        16*(i-1)+1  16*(i-1)+5  16*(i-1)+13  16*(i-1)+17;

        16*(i-1)+1  16*(i-1)+13  16*(i-1)+2  16*(i-1)+7;
        16*(i-1)+2  16*(i-1)+13  16*(i-1)+1  16*(i-1)+5;
        16*(i-1)+5  16*(i-1)+13  16*(i-1)+17  16*(i-1)+18;
        16*(i-1)+7  16*(i-1)+13  16*(i-1)+18  16*(i-1)+17;
        ];
        
end

rotNum=size(rot_spr_4N.node_ijkl_mat);
rotNum=rotNum(1);

rot_spr_4N.rot_spr_K_vec=(10^8)*ones(rotNum,1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;



%% Define 3 Node Rotational Spring
for i=1:N+1    
    rot_spr_3N.node_ijk_mat=[
        rot_spr_3N.node_ijk_mat;
        16*(i-1)+1  16*(i-1)+2  16*(i-1)+4;
        16*(i-1)+2  16*(i-1)+4  16*(i-1)+3;
        16*(i-1)+4  16*(i-1)+3  16*(i-1)+1;
        16*(i-1)+3  16*(i-1)+1  16*(i-1)+2;
        ];
        
end

rot3Num=size(rot_spr_3N.node_ijk_mat,1);
rot_spr_3N.rot_spr_K_vec=(10^8)*ones(rot3Num,1);

plots.Plot_Shape_RotSpr_3N_Number()



%% Initialize Assembly 
assembly.Initialize_Assembly;



%% Calculate self-weight 

% Steel properties
rho_steel=7850;         % Density of steel in kg/m^3
g=9.81;                 % Gravitational acceleration in m/s^2

% Bar elements
A_bar=barA;             % Cross-sectional area of bars in m^2 (matches bar.A_vec)

% Calculate total length of all bars
L_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;

for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    p1=coords(n1,:);
    p2=coords(n2,:);
    len=norm(p1-p2);
    L_total=L_total+len;
end

% Total bar weight (Newtons)
W_bar=A_bar*L_total*rho_steel*g;   


%% Set up solver + Distributed load on bottom nodes (full-length)
nr=Solver_NR_Loading;
nr.assembly=assembly;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

nr.supp=[1 1 1 1;
         2 1 1 1;
         16*N+1 1 1 1;
         16*N+2 1 1 1;];

% force increment of each node per node
force=4000;   % N


for i=1:100

    % Nonlinear solver settings
    nr.increStep=1;
    nr.iterMax=50;
    nr.tol=1e-5; 

    nr.load=[];
    total_F=0;

    for k=1:N-1
        nr.load=[nr.load;
            17+(k-1)*16 0 0 -force*i;
            18+(k-1)*16 0 0 -force*i;
            ];
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