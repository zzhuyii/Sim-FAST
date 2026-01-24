clear all
close all
clc
tic

%% Define Geometry
L=1;
gap=0;
N=8;

node=Elements_Nodes;

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        2*L*(i-1), 0, L;
        2*L*(i-1), 0, 0;
        2*L*(i-1), 2*L, 0;
        2*L*(i-1), 2*L, L;
        2*L*(i-1)+L, 0, L;
        2*L*(i-1)+L, 0, 0;
        2*L*(i-1)+L, L, 0;
        2*L*(i-1)+L, 2*L, 0;
        2*L*(i-1)+L, 2*L, L;];        
end

node.coordinates_mat=[node.coordinates_mat;
        2*L*N, 0, L;
        2*L*N, 0, 0;
        2*L*N, 2*L, 0;
        2*L*N, 2*L, L; ];


%% Define assembly
assembly=Assembly_Origami;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;
bar=Vec_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_4N=rot_spr_4N;

%% Define Plotting Functions
plots=Plot_Kirigami_Truss;
plots.assembly=assembly;
plots.displayRange=[-1; 2*(N+1); -1.5; 3.5; -1; 2];
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
cst.t_vec=0.05*ones(cstNum,1);
cst.E_vec=2*10^9*ones(cstNum,1);
cst.v_vec=0.2*ones(cstNum,1);

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
        ];
end

barNum=size(bar.node_ij_mat,1);
bar.A_vec=0.01*ones(barNum,1);
bar.E_vec=2*10^9*ones(barNum,1);

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
        9*(i-1)+8  9*(i-1)+9  9*(i-1)+12  9*(i-1)+13;        ];
        
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
rot_spr_4N.rot_spr_K_vec=ones(rotNum,1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;

assembly.Initialize_Assembly;


%% Calculate Self-weight of the Bridge
rho_steel = 7850;         % Density of steel in kg/m^3
A_bar = 0.01;             % Cross-sectional area of bar elements in m^2 (consistent with bar.A_vec)

% a. Calculate the total length of all bar elements
L_total = 0;
barNodeMat = bar.node_ij_mat;
coords = node.coordinates_mat;

for i = 1:size(barNodeMat,1)
    n1 = barNodeMat(i,1);
    n2 = barNodeMat(i,2);
    p1 = coords(n1,:);
    p2 = coords(n2,:);
    len = norm(p1 - p2);
    L_total = L_total + len;
end

% b. Calculate the total weight of all bars (in Newtons)
W_bar = A_bar * L_total * rho_steel * 9.81;

% c. Calculate the total area of all CST (triangular plate) elements
A_cst_total = 0;
cstNodeMat = cst.node_ijk_mat;
for i = 1:size(cstNodeMat,1)
    n1 = cstNodeMat(i,1);
    n2 = cstNodeMat(i,2);
    n3 = cstNodeMat(i,3);
    p1 = coords(n1,:);
    p2 = coords(n2,:);
    p3 = coords(n3,:);
    % Use Heron's formula to calculate triangle area
    a = norm(p1 - p2);
    b = norm(p2 - p3);
    c = norm(p3 - p1);
    s = (a + b + c) / 2;
    area = sqrt(s * (s - a) * (s - b) * (s - c));
    A_cst_total = A_cst_total + area;
end

% d. Calculate the total weight of all CST elements
t_cst = 0.005;  % Plate thickness, assumed to be 5 mm (0.005 m)
W_cst = A_cst_total * t_cst * rho_steel * 9.81; % in Newtons

% Total self-weight
W_total = W_bar + W_cst;

fprintf('-----------------------------\n');
fprintf('Total bar length: %.2f m\n', L_total);
fprintf('Total CST area: %.2f m^2\n', A_cst_total);
fprintf('Bar weight: %.2f N\n', W_bar);
fprintf('CST weight: %.2f N\n', W_cst);
fprintf('Total self-weight: %.2f N\n', W_total);
fprintf('-----------------------------\n');


%% Set up solver
sf=Solver_NR_Folding_4N;
sf.assembly=assembly;

for i=1:N+1
    sf.supp=[sf.supp;
         9*(i-1)+2 0 1 1;
         9*(i-1)+3 0 1 1;];
end
sf.supp([1 2],2)=1;

sf.targetRot=pi*ones(size(rot_spr_4N.theta_stress_free_vec));

sf.increStep=100;
sf.iterMax=20;
sf.tol=1*10^-5;

Uhis1=sf.Solve;

sf.increStep=400;
targetFold=0.95*pi;

for i=1:N
    sf.targetRot(16*(i-1)+3)=pi-targetFold;
    sf.targetRot(16*(i-1)+4)=pi-targetFold;
    sf.targetRot(16*(i-1)+7)=pi-targetFold;
    sf.targetRot(16*(i-1)+8)=pi-targetFold;

    sf.targetRot(16*(i-1)+13)=pi+targetFold;
    sf.targetRot(16*(i-1)+14)=pi+targetFold;
    sf.targetRot(16*(i-1)+10)=pi+targetFold;
    sf.targetRot(16*(i-1)+11)=pi+targetFold;
    sf.targetRot(16*(i-1)+15)=pi-targetFold;
    sf.targetRot(16*(i-1)+16)=pi-targetFold;
    sf.targetRot(16*(i-1)+9)=pi+targetFold;
    sf.targetRot(16*(i-1)+12)=pi+targetFold;
end

for i=1:N-1
    sf.targetRot(16*N+3*(i-1)+1)=pi-targetFold;
    sf.targetRot(16*N+3*(i-1)+2)=pi+targetFold;
    sf.targetRot(16*N+3*(i-1)+3)=pi-targetFold;
end

Uhis2=sf.Solve;

toc
plots.Plot_Deformed_Shape(squeeze(Uhis2(end,:,:)))

Uhis=cat(1,Uhis1,Uhis2);

plots.fileName='Origami_Bridge_Deploy.gif';
plots.Plot_Deformed_His(Uhis(1:10:end,:,:))


