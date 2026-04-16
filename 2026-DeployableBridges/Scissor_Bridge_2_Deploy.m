clear all;
clc;
close all;

%% Initialize the scissor 
% Number of Sections
N=8;

% Height of the bridge
H=2; % (m)

% Length of the section
L=2; % (m)

% HSS 8X4X5/16 A500 Grade C Fy=50ksi
barA=0.00415;
barE=2*10^11;

% Second moment of inertia in both direction
Iy=21.2*10^-6; 
Ix=7.16*10^-6; 

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
rotSpr3N.rot_spr_K_vec=kspr*ones(rotNum,1)*10;
rotSpr3N.delta=10^-3;

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


%% Set up solver   
ta=Solver_NR_TrussAction;
ta.assembly = assembly;

ta.supp = [(1:84)'  zeros(84,1)  ones(84,1)  zeros(84,1)];
ta.supp(1,:) =  [ 1   1 1 1];
ta.supp(2,:) =  [ 2   1 1 1];
ta.supp(3,:) =  [ 3   1 1 0];
ta.supp(4,:) =  [ 4   1 1 0];

base_L0=actBar.L0_vec; 
ta.targetL0=base_L0;

for step=1:700

    ta.increStep = 1; 

    dL=0.001*step;
    
    for i=1:actBarNum1
        ta.targetL0(i)=base_L0(i)+dL;
    end
    
    theta=acos((L+dL)/sqrt(2)/L );
    L2=L/sqrt(2)*sin(theta);
    L3 = sqrt((L/2)^2-L2^2);

    for i=(actBarNum1+1):actBarNum
        ta.targetL0(i)=base_L0(i)+dL/2-L3;
    end
    
    ta.iterMax = 40;
    ta.tol = 5*10^-3;
    
    Utemp = ta.Solve();  
    Uhis(step,:,:)=squeeze(Utemp);

    a=1;
end

plots.fileName = 'Scissor_Bridge_2_Deploy.gif';
plots.Plot_Deformed_His(Uhis(1:10:end,:,:));

U_end = squeeze(Uhis(end, :, :));  
plots.Plot_Deformed_Shape(U_end);


