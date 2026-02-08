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

% The three-node rotational spring stiffness
barL=sqrt(H^2+L^2);
kspr=barE*I/barL;



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
        10*(i-1)+13  10*(i-1)+10;

        10*(i-1)+2  10*(i-1)+7;
        10*(i-1)+8  10*(i-1)+11;
        ];
end

i=N+1;
bar.node_ij_mat=[
    bar.node_ij_mat;
    10*(i-1)+1  10*(i-1)+2;
    10*(i-1)+3  10*(i-1)+4;     
    ];

barNum=size(bar.node_ij_mat);
barNum=barNum(1);
bar.A_vec=barA*ones(barNum,1);
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
rotSpr3N.rot_spr_K_vec=kspr*ones(rotNum,1)*100;

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
force=2000;   % N


for i=1:100

    % Nonlinear solver settings
    nr.increStep=1;
    nr.iterMax=50;
    nr.tol=1e-5;  

    nr.load=[];
    total_F=0;
    for k=1:N-1
        nr.load=[nr.load;
            10*(k-1)+1 0 0 -force*i;
            10*(k-1)+2 0 0 -force*i;];
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




