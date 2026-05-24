%% Initialize the solver
clear all;
clc;
close all;


%% Define the Geometry of origami
% Here we generate the geometry of the origami

% Define the nodal coordinate before meshing
R=50*10^(-3);
H=50*10^(-3);
theta=30/180*pi;
N=6;
M=1;

% Stiffness parameters of the structure
sprStiff=0.000001;

barE=10*10^6; % Young's modulus
barH=1*10^6; % Harden modulus
barSigmaY=1*10^6; % Yielding stress


panelv=0.3; % Poisson's ratio for panel
panelt=1*10^-3; % Thickness for panel

% Find the area of bars
alpha=2*pi/N;
node1=[R*cos(theta),R*sin(theta),0];
node2=[R*cos(alpha+theta),R*sin(alpha+theta),0];
node3=[R*cos(alpha+theta*2),R*sin(alpha+theta*2),H];
v1=node2-node1;
v2=node3-node1;
TriangleArea=norm(cross(v2,v1))/2;
TriangleLength=norm(v1)+norm(v2)+norm(node3-node2);

barA=TriangleArea*panelt/TriangleLength*2/(1-panelv);
% The total area of bar is derived through preserving the total volumn
% of panel, which preserve the total elastic energy under uniform dilation
% This formula is from Liu's Merlin2 paper.


% Initialize Elements
bar=Vec_Elements_Bars_Plastic;
rotSpr=Vec_Elements_RotSprings_4N;
node=Elements_Nodes;


%% Geometry of Kresling origami

for i=1:M+1
    for j=1:N
        node.coordinates_mat=[node.coordinates_mat;
            R*cos(j*alpha+i*theta),R*sin(j*alpha+i*theta),(i-1)*H];
    end
end

% Add the bar element that is horizontal
for i=1:M+1
    for j=1:N
        if j ~=N
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i-1)*N+j+1];
        else
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i-1)*N+1];
        end
    end
end

% Add those diagonal bars
for i=1:M
    for j=1:N
        if j ~=N
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i)*N+j];
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i)*N+j+1];
        else
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i)*N+j];
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i)*N+1];
        end
    end
end

% Initialize other parameters
barNum=size(bar.node_ij_mat);
barNum=barNum(1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);
bar.H_vec=barH*ones(barNum,1);
bar.sigma_y_vec=barSigmaY*ones(barNum,1);

bar.alpha_current_vec=zeros(barNum,1);
bar.strain_plastic_current_vec=zeros(barNum,1);


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
assembly=Assembly_Origami();
assembly.node=node;
assembly.bar=bar;
assembly.rotSpr=rotSpr;

assembly.Initialize_Assembly()

%% Plot for investigation
plots=Plot_Origami();
plots.displayRange=0.1;
plots.displayRangeRatio=1;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber()
plots.Plot_Shape_BarNumber()
plots.Plot_Shape_SprNumber()

% set up panels for ploting
panelNum=1;
for i=1:M
    for j=1:N
        if j ~=N
            plots.panelConnection{panelNum}=[
                (i-1)*N+j,(i-1)*N+j+1,(i)*N+j+1];
            panelNum=panelNum+1;
            plots.panelConnection{panelNum}=[
                (i-1)*N+j,(i)*N+j,(i)*N+j+1];
            panelNum=panelNum+1;
        else
            plots.panelConnection{panelNum}=[
                (i-1)*N+j,(i-1)*N+1,(i)*N+1];
            panelNum=panelNum+1;
            plots.panelConnection{panelNum}=[
                (i-1)*N+j,(i)*N+j,(i)*N+1];
            panelNum=panelNum+1;
        end
    end
end



%% Setup the loading controller
increStep=10;

dc=Solver_DC_Plastic;
dc.assembly=assembly;
dc.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;
         4,1,1,1;
         5,1,1,1;
         6,1,1,1];

force=0.1;

dc.selectedRefDisp=[6*M+1,3];

dc.load=[6*M+1,0,0,force;
         6*M+2,0,0,force;
         6*M+3,0,0,force;
         6*M+4,0,0,force;
         6*M+5,0,0,force;
         6*M+6,0,0,force;];

dc.increStep=100;
dc.tol=10^-5;
dc.iterMax=50;

[Uhis,Fhis]=dc.Solve();

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)));

plots.fileName='Kresling.gif';
plots.Plot_DeformedHis(Uhis(1:2:end,:,:));



%% Find the reaction force and loading results

FrefHis=zeros(dc.increStep,1);
UrefHis=zeros(dc.increStep,1);

for i=1:dc.increStep
    UrefHis(i)=Uhis(i,dc.selectedRefDisp(1),dc.selectedRefDisp(2));
    FrefHis(i)=Fhis(i,(dc.selectedRefDisp(1)-1)*3+dc.selectedRefDisp(2));
end

figure
plot([0;UrefHis],[0;FrefHis])
xlabel('Z deformation of top node (m)') 
ylabel('Applied Force (N)') 

