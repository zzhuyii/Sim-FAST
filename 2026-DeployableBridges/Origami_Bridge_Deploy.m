clear all
close all
clc
tic

%% Define Geometry
% Number of Sections
N=4;

% Height of the bridge
L=2;

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
assembly=Assembly_Origami;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;
bar=Vec_Elements_Bars;
node=Elements_Nodes;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_4N=rot_spr_4N;


%% Define Nodal Coordinates
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
        2*L*N, 2*L, L; 
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
plots.Plot_Deformed_Shape(squeeze(Uhis1(end,:,:)))

sf.increStep=2000;
targetFold=0.9*pi;

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

Uhis=cat(1,Uhis1,Uhis2(1:10:end,:,:));

plots.fileName='Origami_Bridge_Deploy.gif';
plots.Plot_Deformed_His(Uhis(1:4:end,:,:))


