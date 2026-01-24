clear all;
clc;
close all;

%% Initialize the scissor 
N=8;
H=1; %0.2286; % (m)
L=1; %0.2286; % (m)
barA=0.0063*0.01; 
barE=2*10^9; % (Pa)

I=1/12*0.01^4; 
barL=sqrt(H^2+L^2);
kspr=3*barE*I/barL;
node=Elements_Nodes();

%% Define the nodal coordinates
for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        L*(i-1), 0, 0;
        L*(i-1), L, 0;
        L*(i-1), 0, L;
        L*(i-1), L, L;

        0.25*L+L*(i-1), 0, 0.25*L;
        0.25*L+L*(i-1), L, 0.25*L;
        0.25*L+L*(i-1), 0, 0.75*L;
        0.25*L+L*(i-1), L, 0.75*L; %8

        0.75*L+L*(i-1), 0, 0.25*L;
        0.75*L+L*(i-1), L, 0.25*L;
        0.75*L+L*(i-1), 0, 0.75*L;
        0.75*L+L*(i-1), L, 0.75*L;

        (i-1)*L+L/2, 0, L/2;
        (i-1)*L+L/2, L, L/2;
        (i-1)*L+L/2, 0, 0;
        (i-1)*L+L/2, L, 0; %16
           
        (i-1)*L+L/2, 0, L;
        (i-1)*L+L/2, L, L;
        (i-1)*L+L/2, L/2, L;
        ];
end

node.coordinates_mat=[node.coordinates_mat;
        L*N, 0, 0;
        L*N, L, 0;
        L*N, 0, L;
        L*N, L, L;
        ];

%% Define assembly
assembly=Assembly_Scissor_Bridge(); 
cst=Vec_Elements_CST;
rotSpr3N=CD_Elements_RotSprings_3N; 
rotSpr4N=Vec_Elements_RotSprings_4N;
bar=Std_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_3N=rotSpr3N;
assembly.rot_spr_4N=rotSpr4N;

%% Define Plotting Functions
plots=Plot_Scissor_Bridge; 
plots.assembly=assembly;
plots.displayRange=[-0.3*L; L*(N+1); -0.3*L; 1.3*L; -0.3*L; 1.5*L]; 
plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_Node_Number;

%% Define Triangle
for i=1:N
cst.node_ijk_mat=[cst.node_ijk_mat;
    19*(i-1)+1  19*(i-1)+2  19*(i-1)+15;
    19*(i-1)+2  19*(i-1)+15 19*(i-1)+16;
    19*(i-1)+16 19*(i-1)+15 19*(i-1)+20;
    19*(i-1)+20 19*(i-1)+16 19*(i-1)+21;
];
end

cstNum=size(cst.node_ijk_mat);
cstNum=cstNum(1);
cst.t_vec=0.0063*ones(cstNum,1); % m
cst.E_vec=2*10^9*ones(cstNum,1); % Pa
cst.v_vec=0.25*ones(cstNum,1); %
plots.Plot_Shape_CST_Number;

%% Define bar
for i=1:N
    bar.node_ij_mat=[
        bar.node_ij_mat;
        19*(i-1)+1   19*(i-1)+5;
        19*(i-1)+5   19*(i-1)+13;
        19*(i-1)+13  19*(i-1)+11;
        19*(i-1)+11  19*(i-1)+22;

        19*(i-1)+3   19*(i-1)+7;
        19*(i-1)+7   19*(i-1)+13;
        19*(i-1)+13  19*(i-1)+9;
        19*(i-1)+9   19*(i-1)+20;

        19*(i-1)+2   19*(i-1)+6;
        19*(i-1)+6   19*(i-1)+14;
        19*(i-1)+14  19*(i-1)+12;
        19*(i-1)+12  19*(i-1)+23;

        19*(i-1)+4   19*(i-1)+8;
        19*(i-1)+8   19*(i-1)+14;
        19*(i-1)+14  19*(i-1)+10;
        19*(i-1)+10  19*(i-1)+21;
        
        % Start of top
        19*(i-1)+3   19*(i-1)+4;
        19*(i-1)+3   19*(i-1)+17;
        19*(i-1)+3   19*(i-1)+19;        
        19*(i-1)+4   19*(i-1)+19;
        19*(i-1)+4   19*(i-1)+18;

        19*(i-1)+17  19*(i-1)+19;
        19*(i-1)+18  19*(i-1)+19;
        
        19*(i-1)+18  19*(i-1)+23;               
        19*(i-1)+19   19*(i-1)+23;
        19*(i-1)+17  19*(i-1)+22; 
        19*(i-1)+19  19*(i-1)+22;

        % add new bars
        19*(i-1)+1    19*(i-1)+2;
        19*(i-1)+15   19*(i-1)+16;
        19*(i-1)+1    19*(i-1)+15;
        19*(i-1)+15   19*(i-1)+20;
        19*(i-1)+2    19*(i-1)+16;
        19*(i-1)+16   19*(i-1)+21;
        
        19*(i-1)+2    19*(i-1)+15;
        19*(i-1)+16   19*(i-1)+20;

        ];
end

bar.node_ij_mat=[
        bar.node_ij_mat;
        19*N+3   19*N+4;
        19*N+1   19*N+2;
        ];

barNum=size(bar.node_ij_mat);
barNum=barNum(1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);
plots.Plot_Shape_Bar_Number();

%% Define 3 Node Rotational Spring
for i=1:N    
    rotSpr3N.node_ijk_mat = [rotSpr3N.node_ijk_mat;
         19*(i-1)+1  19*(i-1)+5  19*(i-1)+13;
         19*(i-1)+5  19*(i-1)+13 19*(i-1)+11;
         19*(i-1)+13 19*(i-1)+11 19*(i-1)+22;

         19*(i-1)+3  19*(i-1)+7  19*(i-1)+13;
         19*(i-1)+7  19*(i-1)+13 19*(i-1)+9;
         19*(i-1)+13 19*(i-1)+9  19*(i-1)+20;

         19*(i-1)+2  19*(i-1)+6  19*(i-1)+14;
         19*(i-1)+6  19*(i-1)+14 19*(i-1)+12;
         19*(i-1)+14 19*(i-1)+12 19*(i-1)+23;

         19*(i-1)+4  19*(i-1)+8  19*(i-1)+14;
         19*(i-1)+8  19*(i-1)+14 19*(i-1)+10;
         19*(i-1)+14 19*(i-1)+10 19*(i-1)+21;
         ];
end

rotNum=size(rotSpr3N.node_ijk_mat,1);
rotSpr3N.rot_spr_K_vec=kspr*ones(rotNum,1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_RotSpr_3N_Number;



%% Set up four node rotational spring
for i=1:N    
    rotSpr4N.node_ijkl_mat = [rotSpr4N.node_ijkl_mat;
         19*(i-1)+1   19*(i-1)+2   19*(i-1)+15  19*(i-1)+16;
         19*(i-1)+2   19*(i-1)+15  19*(i-1)+16  19*(i-1)+20;
         19*(i-1)+15  19*(i-1)+16  19*(i-1)+20  19*(i-1)+21;

         19*(i-1)+3  19*(i-1)+17  19*(i-1)+19  19*(i-1)+22;
         19*(i-1)+4  19*(i-1)+18  19*(i-1)+19  19*(i-1)+23;

         19*(i-1)+3  19*(i-1)+4  19*(i-1)+19  19*(i-1)+18;
         19*(i-1)+4  19*(i-1)+3  19*(i-1)+19  19*(i-1)+17;
         19*(i-1)+18  19*(i-1)+19  19*(i-1)+23  19*(i-1)+22;
         19*(i-1)+17  19*(i-1)+22  19*(i-1)+19  19*(i-1)+23;
         ];
end

rotNum4N=size(rotSpr4N.node_ijkl_mat,1);
rotSpr4N.rot_spr_K_vec=0.01*ones(rotNum4N,1);
plots.Plot_Shape_RotSpr_4N_Number;
assembly.Initialize_Assembly;

%% Set up solver   
sf = Solver_NR_Folding_4N;
sf.assembly = assembly;

sf.supp = [1    1 1 1;
           2    1 1 1;
           3    1 1 0;
           4    1 1 0;
           ];

sf.increStep = 50; %200;
sf.targetRot = assembly.rot_spr_4N.theta_current_vec;

targetAngle=0.9*pi;

for i=1:N
    sf.targetRot(9*(i-1)+2)=pi-targetAngle;
    sf.targetRot(9*(i-1)+4)=pi+targetAngle;
    sf.targetRot(9*(i-1)+5)=pi+targetAngle;
end

sf.iterMax = 20;
sf.tol = 10^-5;

Uhis = sf.Solve();  
plots.fileName = 'Scissor_Bridge_Deploy.gif';
plots.Plot_Deformed_His(Uhis(1:10:end,:,:));

U_end = squeeze(Uhis(end, :, :));  
plots.Plot_Deformed_Shape(U_end);