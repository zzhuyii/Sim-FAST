clear all;
clc;
close all;

cd('G:\My Drive\Dr. Zhu\Code\Sim-FAST\2026-DeployableBridges');

%% Initialize the scissor 
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
        8*(i-1)+1  8*(i-1)+2  8*(i-1)+7;
        8*(i-1)+2  8*(i-1)+7 8*(i-1)+8;
        8*(i-1)+7 8*(i-1)+8 8*(i-1)+9;
        8*(i-1)+8 8*(i-1)+10 8*(i-1)+9;
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
        8*(i-1)+1   8*(i-1)+7;
        8*(i-1)+7   8*(i-1)+9;
        8*(i-1)+2   8*(i-1)+8;
        8*(i-1)+8   8*(i-1)+10;

        8*(i-1)+3   8*(i-1)+4;
        8*(i-1)+9   8*(i-1)+10;
        8*(i-1)+1   8*(i-1)+2;
        8*(i-1)+7   8*(i-1)+8;

        8*(i-1)+1  8*(i-1)+5;
        8*(i-1)+3  8*(i-1)+5;
        8*(i-1)+2  8*(i-1)+6;
        8*(i-1)+4  8*(i-1)+6;  

        8*(i-1)+5  8*(i-1)+9;
        8*(i-1)+5  8*(i-1)+11;
        8*(i-1)+6  8*(i-1)+10;
        8*(i-1)+6  8*(i-1)+12;

        8*(i-1)+2  8*(i-1)+7;
        8*(i-1)+8  8*(i-1)+9;
        ];
end

i=N+1;
bar.node_ij_mat=[
    bar.node_ij_mat;
    8*(i-1)+1  8*(i-1)+2;
    8*(i-1)+3  8*(i-1)+4;        
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
         8*(i-1)+1  8*(i-1)+5  8*(i-1)+11;
         8*(i-1)+3  8*(i-1)+5  8*(i-1)+9;
         8*(i-1)+2  8*(i-1)+6  8*(i-1)+12;
         8*(i-1)+4  8*(i-1)+6  8*(i-1)+10;

         8*(i-1)+4  8*(i-1)+3  8*(i-1)+5;
         8*(i-1)+3  8*(i-1)+4  8*(i-1)+6;
         8*(i-1)+2  8*(i-1)+1  8*(i-1)+5;
         8*(i-1)+6  8*(i-1)+2  8*(i-1)+1;

         8*(i-1)+5  8*(i-1)+11  8*(i-1)+12;
         8*(i-1)+11  8*(i-1)+12  8*(i-1)+6;
         8*(i-1)+5  8*(i-1)+9  8*(i-1)+10;
         8*(i-1)+9  8*(i-1)+10  8*(i-1)+6;

         ];
end

rotNum=size(rotSpr3N.node_ijk_mat,1);
rotSpr3N.rot_spr_K_vec=kspr*ones(rotNum,1)*100;
rotSpr3N.delta=10^-3;

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_RotSpr_3N_Number;



%% Set up four node rotational spring
for i=1:N    
    rotSpr4N.node_ijkl_mat = [rotSpr4N.node_ijkl_mat;
         8*(i-1)+1   8*(i-1)+2   8*(i-1)+7  8*(i-1)+8;
         8*(i-1)+7  8*(i-1)+8  8*(i-1)+9  8*(i-1)+10;
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
        8*(i-1)+1 8*(i-1)+3;
        8*(i-1)+2 8*(i-1)+4;
        ];
end

i=N+1;
actBar.node_ij_mat=[
    actBar.node_ij_mat;
    8*(i-1)+1 8*(i-1)+3;
    8*(i-1)+2 8*(i-1)+4;
    ];

actBarNum1=size(actBar.node_ij_mat,1);

for i=1:N
    actBar.node_ij_mat=[
        actBar.node_ij_mat;
        8*(i-1)+7 8*(i-1)+5;
        8*(i-1)+8 8*(i-1)+6;
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

ta.supp = [(1:68)'  zeros(68,1)  ones(68,1)  zeros(68,1)];
ta.supp(1,:) =  [ 1   1 1 1];
ta.supp(2,:) =  [ 2   1 1 1];
ta.supp(3,:) =  [ 3   1 1 0];
ta.supp(4,:) =  [ 4   1 1 0];
ta.supp(9,:) =  [ 9   0 1 1];
ta.supp(10,:) =  [ 10   0 1 1];
ta.supp(17,:) =  [ 17   0 1 1];
ta.supp(18,:) =  [ 18   0 1 1];
ta.supp(25,:) =  [ 25   0 1 1];
ta.supp(26,:) =  [ 26   0 1 1];
ta.supp(33,:) =  [ 33   0 1 1];
ta.supp(34,:) =  [ 34   0 1 1];
ta.supp(41,:) =  [ 41   0 1 1];
ta.supp(42,:) =  [ 42   0 1 1];
ta.supp(49,:) =  [ 49   0 1 1];
ta.supp(50,:) =  [ 50   0 1 1];
ta.supp(65,:) =  [ 65   0 1 1];
ta.supp(66,:) =  [ 66   0 1 1];



% ta.supp = [1    1 1 1;
%            2    1 1 1;
%            3    1 1 0;
%            4    1 1 0;
%            65   0 0 1;
%            66   0 0 1;
%            ];

base_L0=actBar.L0_vec; 
ta.targetL0=base_L0;

nodeNum = size(node.coordinates_mat, 1);
Uhis = zeros(2000, nodeNum, 3);

for step=1:2000

    ta.increStep = 1; 

    % dL=0.001*step; 

    if step<=200
        dL=0.001*step; 
    elseif step <= 400
        dL=0.001*200+0.0004*(step-200);
    else
        dL=0.001*200+0.0004*(400-200)+0.0001*(step-400);
    end

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
    ta.tol = 10^-4;

    Utemp = ta.Solve();  
    Uhis(step,:,:)=squeeze(Utemp);

    a=1;
end


plots.fileName = 'Scissor_Bridge_Deploy.gif';
plots.Plot_Deformed_His(Uhis(1:20:end,:,:));
U_end = squeeze(Uhis(end, :, :));  
plots.Plot_Deformed_Shape(U_end);

% Store the deformation history
save('ScissorUhis.mat', 'Uhis');
UhisNew=load('ScissorUhis.mat');
UhisNew=UhisNew.Uhis;


% 
% 
% plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)))
% plots.fileName = 'Scissor_Bridge_Deploy.gif';
% plots.Plot_Deformed_His(Uhis(1:20:end,:,:));
% 
% % store the deformation history
% save('ScissorUhis.mat','Uhis'); % Saves to a .mat file
% U_end = squeeze(Uhis(end, :, :));  
% plots.Plot_Deformed_Shape(U_end);
% 
% base_L0=actBar.L0_vec; 
% ta.targetL0=base_L0;
% 
% % dL_final: total actuator elongation at full deployment (same as 2000-step schedule)
% dL_final = 0.001*200 + 0.0004*200 + 0.0001*1600;  % = 0.44 m
% 
% for step=1:300
% 
%     ta.increStep = 1;
% 
%     % Linear ramp across 10 steps to reach the same final dL
%     dL = dL_final / 300 * step;
% 
%     for i=1:actBarNum1
%         ta.targetL0(i)=base_L0(i)+dL;
%     end
% 
%     theta=acos((L+dL)/sqrt(2)/L );
%     L2=L/sqrt(2)*sin(theta);
%     L3 = sqrt((L/2)^2-L2^2);
% 
%     for i=(actBarNum1+1):actBarNum
%         ta.targetL0(i)=base_L0(i)+dL/2-L3;
%     end
% 
%     ta.iterMax = 40;
%     ta.tol = 10^-4;
% 
%     Utemp = ta.Solve();  
%     Uhis(step,:,:)=squeeze(Utemp);
% 
%     a=1;
% end