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
% HSS 8X4X5/16 A500 Grade C Fy=50ksi
barA=0.00415;
barE=2*10^11;

% Second moment of inertia in both direction
Iy=21.2*10^-6; 
Ix=7.16*10^-6;  

% We assume a soft panel so that only the truss is taking global load
% Thus, panel Young's modulus is 20 MPa
panel_E=2*10^7;
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

rot_spr_4N.rot_spr_K_vec=100*ones(rotNum,1);

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



%% Set up solver
sf=Solver_NR_Folding_4N;
sf.assembly=assembly;

sf.supp=[1 1 1 1;
         2 1 1 1;
         3 1 1 1;
         4 1 1 1;];

sf.targetRot=rot_spr_4N.theta_stress_free_vec;

sf.increStep=1000;
sf.iterMax=30;
sf.tol=1*10^-4;

rate=0.9;

for i=1:N
    sf.targetRot((i-1)*24+1)=pi+rate*pi;
    sf.targetRot((i-1)*24+2)=pi+rate*pi;

    sf.targetRot((i-1)*24+3)=pi-rate*pi;
    sf.targetRot((i-1)*24+4)=pi-rate*pi;
    sf.targetRot((i-1)*24+5)=pi-rate*pi;
    sf.targetRot((i-1)*24+6)=pi-rate*pi;

    sf.targetRot((i-1)*24+7)=pi-rate*pi;
    sf.targetRot((i-1)*24+8)=pi-rate*pi;

    sf.targetRot((i-1)*24+9)=pi+rate*pi;
    sf.targetRot((i-1)*24+10)=pi+rate*pi;
    sf.targetRot((i-1)*24+11)=pi+rate*pi;
    sf.targetRot((i-1)*24+12)=pi+rate*pi;

    sf.targetRot((i-1)*24+13)=pi-rate*pi;
    sf.targetRot((i-1)*24+14)=pi-rate*pi;

    sf.targetRot((i-1)*24+15)=pi+rate*pi;
    sf.targetRot((i-1)*24+16)=pi+rate*pi;
    sf.targetRot((i-1)*24+17)=pi+rate*pi;
    sf.targetRot((i-1)*24+18)=pi+rate*pi;

    sf.targetRot((i-1)*24+19)=pi+rate*pi;
    sf.targetRot((i-1)*24+20)=pi+rate*pi;

    sf.targetRot((i-1)*24+21)=pi-rate*pi;
    sf.targetRot((i-1)*24+22)=pi-rate*pi;
    sf.targetRot((i-1)*24+23)=pi-rate*pi;
    sf.targetRot((i-1)*24+24)=pi-rate*pi;
    
end

Uhis=sf.Solve;

toc
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)))
plots.fileName='Kirigami_Truss_Deploy.gif';
plots.Plot_Deformed_His(Uhis(1:20:end,:,:))


% store the deformation history
save('KirigamiUhis.mat','Uhis'); % Saves to a .mat file
UhisNew=load('KirigamiUhis.mat');
UhisNew=UhisNew.Uhis;