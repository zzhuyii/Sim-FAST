%% Initialize the solver
clear all;
clc;
close all;

tic

%% Define the Geometry of origami
% Here we generate the geometry of the origami

% Define the nodal coordinate before meshing
R=50*10^(-3);
H=50*10^(-3);
theta=30/180*pi;
N=6;
M=1;

% Stiffness parameters of the structure
sprStiff=0.00001;

panelE=10*10^6;
panelv=0.3;
panelt=1*10^-3;

cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings_4N;
node=Elements_Nodes;



%% Geometry of Kresling origami

alpha=2*pi/N;
for i=1:M+1
    for j=1:N
        node.coordinates_mat=[node.coordinates_mat;
            R*cos(j*alpha+i*theta),R*sin(j*alpha+i*theta),(i-1)*H];
    end
end


% set up panels for plotting
panelNum=0;
for i=1:M
    for j=1:N
        if j ~=N
            cst.node_ijk_mat=[cst.node_ijk_mat;
                (i-1)*N+j,(i-1)*N+j+1,(i)*N+j+1];
            panelNum=panelNum+1;
            cst.node_ijk_mat=[cst.node_ijk_mat;
                (i-1)*N+j,(i)*N+j,(i)*N+j+1];
            panelNum=panelNum+1;
        else
            cst.node_ijk_mat=[cst.node_ijk_mat;
                (i-1)*N+j,(i-1)*N+1,(i)*N+1];
            panelNum=panelNum+1;
            cst.node_ijk_mat=[cst.node_ijk_mat;
                (i-1)*N+j,(i)*N+j,(i)*N+1];
            panelNum=panelNum+1;
        end
    end
end

cst.v_vec=panelv*ones(panelNum,1);
cst.E_vec=panelE*ones(panelNum,1);
cst.t_vec=panelt*ones(panelNum,1);


%% Set up the rotational springs
% Diagonal rotational springs
for i=1:M
    for j=1:N
        if j ==1
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+N,(i-1)*N+j,(i)*N+j,(i)*N+j+1];
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i)*N+j,(i-1)*N+j,(i)*N+j+1,(i-1)*N+j+1];
        elseif j~=N
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+j-1,(i-1)*N+j,(i)*N+j,(i)*N+j+1];
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i)*N+j,(i-1)*N+j,(i)*N+j+1,(i-1)*N+j+1];
        else
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+j-1,(i-1)*N+j,(i)*N+j,(i)*N+1];
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i)*N+j,(i-1)*N+j,(i)*N+1,(i-1)*N+1];
        end
    end
end

for i=1:M-1
    for j=1:N
        if j~=N
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+j,i*N+j,i*N+j+1,(i+1)*N+j+1];
        else
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+j,i*N+j,i*N+1,(i+1)*N+1];
        end
    end
end

rotSpr.rot_spr_K_vec=sprStiff*ones(M*(2*N)+N*(M-1),1);



%% Initialize assembly
assembly=Assembly_CST_Origami();
assembly.node=node;
assembly.cst=cst;
assembly.rotSpr=rotSpr;

assembly.Initialize_Assembly()

%% Plot for investigation
plots=Plot_CST_Origami();
plots.displayRange=0.1;
plots.displayRangeRatio=1;
plots.assembly=assembly;


plots.Plot_Shape_CSTNumber()
plots.Plot_Shape_NodeNumber()
plots.Plot_Shape_SprNumber()



%% Setup the loading controller
dc=Solver_DC;
dc.assembly=assembly;
dc.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;
         4,1,1,1;
         5,1,1,1;
         6,1,1,1];

force=-1;

dc.selectedRefDisp=[6*M+1,3];

dc.load=[6*M+1,0,0,force;
         6*M+2,0,0,force;
         6*M+3,0,0,force;
         6*M+4,0,0,force;
         6*M+5,0,0,force;
         6*M+6,0,0,force;];

dc.increStep=20;
dc.tol=10^-6;
dc.iterMax=50;

Uhis=dc.Solve();

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)));

% plots.fileName='CST_Kresling.gif';
% plots.Plot_DeformedHis(Uhis);

%% Find the reaction force and loading results

forceHis=zeros(dc.increStep,1);
UrefHis=zeros(dc.increStep,1);

for i=1:dc.increStep
    [F,K]=assembly.Solve_FK(squeeze(Uhis(i,:,:)));
    UrefHis(i)=Uhis(i,dc.selectedRefDisp(1),dc.selectedRefDisp(2));
    forceHis(i)=F(dc.selectedRefDisp(2)+(dc.selectedRefDisp(1)-1)*3);
end

figure
UrefHis=[0;UrefHis];
forceHis=[0;forceHis];
plot(-UrefHis,forceHis)
xlabel('Z deformation of top node (m)') 
ylabel('Applied Force (N)') 

toc

