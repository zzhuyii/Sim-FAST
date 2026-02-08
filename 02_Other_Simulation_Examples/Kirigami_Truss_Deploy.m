clear all
close all
clc
tic

%% Define Geometry
L=1;
w=0.1;
gap=0;
N=4;

node=Elements_Nodes;
node.coordinates_mat=[node.coordinates_mat;
        -w, 0, 0;
        -w, L, 0;
        -w, 0, L;
        -w, L, L;];

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        (w+L)*(i-1), 0, 0;
        (w+L)*(i-1), L, 0;
        (w+L)*(i-1), 0, L;
        (w+L)*(i-1), L, L;

        (w+L)*(i-1)+L/2, 0, 0;
        (w+L)*(i-1)+L/2, 0, gap;
        (w+L)*(i-1)+L/2, L, 0;
        (w+L)*(i-1)+L/2, L, gap;

        (w+L)*(i-1)+L/2, 0, L;
        (w+L)*(i-1)+L/2, gap, L;
        (w+L)*(i-1)+L/2, L, L;
        (w+L)*(i-1)+L/2, L, L-gap;

        (w+L)*(i-1)+L/2, L/2, 0;
        (w+L)*(i-1)+L/2, L/2, L;
        (w+L)*(i-1)+L/2, 0, L/2;
        (w+L)*(i-1)+L/2, L, L/2;

        (w+L)*(i-1)+L, 0, 0;
        (w+L)*(i-1)+L, L, 0;
        (w+L)*(i-1)+L, 0, L;
        (w+L)*(i-1)+L, L, L;];
end

node.coordinates_mat=[node.coordinates_mat;
        (w+L)*N, 0, 0;
        (w+L)*N, L, 0;
        (w+L)*N, 0, L;
        (w+L)*N, L, L;];


%% Define assembly
assembly=Assembly_Kirigami_Truss;
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
plots.displayRange=[-1; L*(N+1); -1; 2; -1; 2];
plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_Node_Number;


%% Define Triangle
for i=1:N+1
     cst.node_ijk_mat=[cst.node_ijk_mat;
        20*(i-1)+1    20*(i-1)+3    20*(i-1)+7;
        20*(i-1)+1    20*(i-1)+5    20*(i-1)+7;
        20*(i-1)+3    20*(i-1)+4    20*(i-1)+8;
        20*(i-1)+3    20*(i-1)+7    20*(i-1)+8;
        20*(i-1)+2    20*(i-1)+4    20*(i-1)+6;
        20*(i-1)+4    20*(i-1)+6    20*(i-1)+8;
        20*(i-1)+3    20*(i-1)+4    20*(i-1)+8;
        20*(i-1)+1    20*(i-1)+2    20*(i-1)+6;
        20*(i-1)+1    20*(i-1)+5    20*(i-1)+6;];
end

for i=1:N
     cst.node_ijk_mat=[cst.node_ijk_mat;
        20*(i-1)+5    20*(i-1)+6    20*(i-1)+17;
        20*(i-1)+5    20*(i-1)+9    20*(i-1)+17;
        20*(i-1)+11    20*(i-1)+6    20*(i-1)+17;
        20*(i-1)+9    20*(i-1)+21    20*(i-1)+17;
        20*(i-1)+11    20*(i-1)+22    20*(i-1)+17;
        20*(i-1)+21    20*(i-1)+22    20*(i-1)+17;];
end

cstNum=size(cst.node_ijk_mat,1);
cst.t_vec=0.05*ones(cstNum,1);
cst.E_vec=2*10^9*ones(cstNum,1);
cst.v_vec=0.2*ones(cstNum,1);

plots.Plot_Shape_CST_Number;

%% Define bar
for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        20*(i-1)+5   20*(i-1)+10;
        20*(i-1)+19  20*(i-1)+10;
        20*(i-1)+5   20*(i-1)+19;
        20*(i-1)+7   20*(i-1)+13;
        20*(i-1)+7   20*(i-1)+19;
        20*(i-1)+13   20*(i-1)+19;
        20*(i-1)+13   20*(i-1)+23;
        20*(i-1)+19   20*(i-1)+23;
        20*(i-1)+19   20*(i-1)+21;
        20*(i-1)+10   20*(i-1)+21; %10

        20*(i-1)+7   20*(i-1)+14;
        20*(i-1)+7   20*(i-1)+18;
        20*(i-1)+14   20*(i-1)+18;
        20*(i-1)+8   20*(i-1)+15;
        20*(i-1)+15   20*(i-1)+18;
        20*(i-1)+8   20*(i-1)+18;
        20*(i-1)+14   20*(i-1)+23;
        20*(i-1)+18   20*(i-1)+23;
        20*(i-1)+18   20*(i-1)+24;
        20*(i-1)+15   20*(i-1)+24; %20

        20*(i-1)+8   20*(i-1)+16;
        20*(i-1)+16   20*(i-1)+20;
        20*(i-1)+8   20*(i-1)+20;
        20*(i-1)+24   20*(i-1)+16;
        20*(i-1)+24   20*(i-1)+20;
        20*(i-1)+6   20*(i-1)+20;
        20*(i-1)+6   20*(i-1)+12;
        20*(i-1)+12   20*(i-1)+20;
        20*(i-1)+22   20*(i-1)+12;
        20*(i-1)+22   20*(i-1)+20; %30

        20*(i-1)+5   20*(i-1)+9;
        20*(i-1)+9   20*(i-1)+17;
        20*(i-1)+5   20*(i-1)+17;
        20*(i-1)+6   20*(i-1)+11;
        20*(i-1)+11   20*(i-1)+17;
        20*(i-1)+6   20*(i-1)+17;
        20*(i-1)+11   20*(i-1)+22;
        20*(i-1)+17   20*(i-1)+22;
        20*(i-1)+9   20*(i-1)+21;
        20*(i-1)+17   20*(i-1)+21; %40
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
        20*(i-1)+5  20*(i-1)+1  20*(i-1)+7  20*(i-1)+2;
        20*(i-1)+1  20*(i-1)+7  20*(i-1)+3  20*(i-1)+8;
        20*(i-1)+7  20*(i-1)+3  20*(i-1)+8  20*(i-1)+4;
        20*(i-1)+3  20*(i-1)+4  20*(i-1)+8  20*(i-1)+6;
        20*(i-1)+2  20*(i-1)+4  20*(i-1)+6  20*(i-1)+8;
        20*(i-1)+4  20*(i-1)+2  20*(i-1)+6  20*(i-1)+1;
        20*(i-1)+2  20*(i-1)+6  20*(i-1)+1  20*(i-1)+5;
        20*(i-1)+6  20*(i-1)+1  20*(i-1)+5  20*(i-1)+7; %8

        20*(i-1)+1  20*(i-1)+7  20*(i-1)+5  20*(i-1)+19;
        20*(i-1)+5  20*(i-1)+7  20*(i-1)+19  20*(i-1)+13;
        20*(i-1)+7  20*(i-1)+13  20*(i-1)+19  20*(i-1)+23;
        20*(i-1)+19  20*(i-1)+23  20*(i-1)+21  20*(i-1)+27;
        20*(i-1)+5  20*(i-1)+10  20*(i-1)+19  20*(i-1)+21;
        20*(i-1)+7  20*(i-1)+5  20*(i-1)+19  20*(i-1)+10;
        20*(i-1)+13  20*(i-1)+19  20*(i-1)+23  20*(i-1)+21;
        20*(i-1)+10  20*(i-1)+19  20*(i-1)+21  20*(i-1)+23; %16

        20*(i-1)+3  20*(i-1)+8  20*(i-1)+7  20*(i-1)+18;
        20*(i-1)+8  20*(i-1)+7  20*(i-1)+18  20*(i-1)+14;
        20*(i-1)+7  20*(i-1)+8  20*(i-1)+18  20*(i-1)+15;
        20*(i-1)+8  20*(i-1)+15  20*(i-1)+18  20*(i-1)+24;
        20*(i-1)+7  20*(i-1)+14  20*(i-1)+18  20*(i-1)+23;
        20*(i-1)+14  20*(i-1)+18  20*(i-1)+23  20*(i-1)+24;
        20*(i-1)+15  20*(i-1)+18  20*(i-1)+24  20*(i-1)+23;
        20*(i-1)+18  20*(i-1)+23  20*(i-1)+24  20*(i-1)+28; %24

        20*(i-1)+4  20*(i-1)+6  20*(i-1)+8  20*(i-1)+20;
        20*(i-1)+6  20*(i-1)+8  20*(i-1)+20  20*(i-1)+16;
        20*(i-1)+8  20*(i-1)+6  20*(i-1)+20  20*(i-1)+12;
        20*(i-1)+8  20*(i-1)+16  20*(i-1)+20  20*(i-1)+24;
        20*(i-1)+6  20*(i-1)+12  20*(i-1)+20  20*(i-1)+22;
        20*(i-1)+24  20*(i-1)+20  20*(i-1)+22  20*(i-1)+12;
        20*(i-1)+16  20*(i-1)+20  20*(i-1)+24  20*(i-1)+22;
        20*(i-1)+26  20*(i-1)+24  20*(i-1)+22  20*(i-1)+20; %32

        20*(i-1)+1  20*(i-1)+5  20*(i-1)+6  20*(i-1)+17;
        20*(i-1)+6  20*(i-1)+5  20*(i-1)+17  20*(i-1)+9;
        20*(i-1)+5  20*(i-1)+6  20*(i-1)+17  20*(i-1)+11;
        20*(i-1)+5  20*(i-1)+9  20*(i-1)+17  20*(i-1)+21;
        20*(i-1)+6  20*(i-1)+17  20*(i-1)+11  20*(i-1)+22;
        20*(i-1)+11  20*(i-1)+17  20*(i-1)+22  20*(i-1)+21;
        20*(i-1)+9  20*(i-1)+17  20*(i-1)+21  20*(i-1)+22;
        20*(i-1)+17  20*(i-1)+21  20*(i-1)+22  20*(i-1)+26; %40
        ];
        
end

rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    20*N+5  20*N+1  20*N+7  20*N+2;
    20*N+1  20*N+7  20*N+3  20*N+8;
    20*N+7  20*N+3  20*N+8  20*N+4;
    20*N+3  20*N+4  20*N+8  20*N+6;
    20*N+2  20*N+4  20*N+6  20*N+8;
    20*N+4  20*N+2  20*N+6  20*N+1;
    20*N+2  20*N+6  20*N+1  20*N+5;
    20*N+6  20*N+1  20*N+5  20*N+7;];


rotNum=size(rot_spr_4N.node_ijkl_mat);
rotNum=rotNum(1);

rot_spr_4N.rot_spr_K_vec=ones(rotNum,1);

factor=100;
for i=1:N+1
    rot_spr_4N.rot_spr_K_vec((i-1)*40+2)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+2);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+4)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+4);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+6)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+6);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+8)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+8);
end

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;

assembly.Initialize_Assembly;


%% Set up solver
sf=Solver_NR_Folding_4N;
sf.assembly=assembly;

sf.supp=[1 1 1 1;
         2 1 1 1;
         3 1 1 1;
         4 1 1 1;];

sf.targetRot=rot_spr_4N.theta_stress_free_vec;

sf.increStep=500;
sf.iterMax=20;
sf.tol=1*10^-4;

rate=0.9;

for i=1:N
    sf.targetRot((i-1)*40+11)=pi+rate*pi;
    sf.targetRot((i-1)*40+13)=pi+rate*pi;

    sf.targetRot((i-1)*40+10)=pi-rate*pi;
    sf.targetRot((i-1)*40+14)=pi-rate*pi;
    sf.targetRot((i-1)*40+15)=pi-rate*pi;
    sf.targetRot((i-1)*40+16)=pi-rate*pi;

    sf.targetRot((i-1)*40+20)=pi-rate*pi;
    sf.targetRot((i-1)*40+21)=pi-rate*pi;

    sf.targetRot((i-1)*40+18)=pi+rate*pi;
    sf.targetRot((i-1)*40+19)=pi+rate*pi;
    sf.targetRot((i-1)*40+22)=pi+rate*pi;
    sf.targetRot((i-1)*40+23)=pi+rate*pi;

    sf.targetRot((i-1)*40+36)=pi+rate*pi;
    sf.targetRot((i-1)*40+37)=pi+rate*pi;

    sf.targetRot((i-1)*40+34)=pi-rate*pi;
    sf.targetRot((i-1)*40+35)=pi-rate*pi;
    sf.targetRot((i-1)*40+38)=pi-rate*pi;
    sf.targetRot((i-1)*40+39)=pi-rate*pi;

    sf.targetRot((i-1)*40+28)=pi-rate*pi;
    sf.targetRot((i-1)*40+29)=pi-rate*pi;

    sf.targetRot((i-1)*40+26)=pi+rate*pi;
    sf.targetRot((i-1)*40+27)=pi+rate*pi;
    sf.targetRot((i-1)*40+30)=pi+rate*pi;
    sf.targetRot((i-1)*40+31)=pi+rate*pi;
    
end

Uhis=sf.Solve;

toc
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)))
plots.fileName='Kirigami_Truss_Deploy.gif';
plots.Plot_Deformed_His(Uhis(1:5:end,:,:))


