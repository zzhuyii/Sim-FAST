%% Initialize the solver
clear all;
clc;
close all;


%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% Define the nodal coordinate before meshing
a=20*10^(-3);
L=50*10^(-3);
lambda=30/180*pi;
omega=60/180*pi;

% Stiffness parameters of the structure
sprStiff=0.000001;
stiffFactor=1000;
barA=10*10^(-3)*10*10^(-3); 
barE=2*10^9; % Young's modulus

% perturbation Level 
perturbationLevel=0*10^(-3);

data=zeros(30,10);

bar=Vec_Elements_Bars;
rotSpr=Vec_Elements_RotSprings_4N;
node=Elements_Nodes;

bar.A_vec=zeros(5,1);
bar.E_vec=zeros(5,1);
bar.L0_vec=zeros(5,1);


%% Bennett Linkage Geometry
% First Links
node.coordinates_mat(1,:)=[L*cos(omega),0,0];
node.coordinates_mat(2,:)=[L*cos(omega)+a*cos(lambda)/sin(omega),0,a*sin(lambda)];
node.coordinates_mat(3,:)=[L*cos(omega)+a*cos(lambda)/sin(omega)+a*sin(lambda)/sin(omega),0,a*sin(lambda)-a*cos(lambda)];
node.coordinates_mat(4,:)=[L*cos(omega)+a*sin(lambda)/sin(omega),0,-a*cos(lambda)];

bar.E_vec(1:5)=barE;
bar.A_vec(1:5)=barA;
bar.node_ij_mat(1,:)=[1,2];
bar.node_ij_mat(2,:)=[2,3];
bar.node_ij_mat(3,:)=[3,4];
bar.node_ij_mat(4,:)=[1,4];
bar.node_ij_mat(5,:)=[1,3];

node.coordinates_mat(5,:)=[0,L*sin(omega),0];
node.coordinates_mat(6,:)=[0,L*sin(omega)+a*cos(lambda)/cos(omega),a*sin(lambda)];
node.coordinates_mat(7,:)=[0,L*sin(omega)+a*cos(lambda)/cos(omega)+a*sin(lambda)/cos(omega),a*sin(lambda)-a*cos(lambda)];
node.coordinates_mat(8,:)=[0,L*sin(omega)+a*sin(lambda)/cos(omega),-a*cos(lambda)];

bar.E_vec(6:10)=barE;
bar.A_vec(6:10)=barA;
bar.node_ij_mat(6,:)=[5,6];
bar.node_ij_mat(7,:)=[6,7];
bar.node_ij_mat(8,:)=[7,8];
bar.node_ij_mat(9,:)=[5,8];
bar.node_ij_mat(10,:)=[5,7];

bar.E_vec(11:14)=barE;
bar.A_vec(11:14)=barA;
bar.node_ij_mat(11,:)=[1,5];
bar.node_ij_mat(12,:)=[2,6];
bar.node_ij_mat(13,:)=[3,7];
bar.node_ij_mat(14,:)=[4,8];

bar.E_vec(15:18)=barE;
bar.A_vec(15:18)=barA;
bar.node_ij_mat(15,:)=[1,6];
bar.node_ij_mat(16,:)=[2,7];
bar.node_ij_mat(17,:)=[3,8];
bar.node_ij_mat(18,:)=[4,5];

% Second Link
node.coordinates_mat(9,:)=[0,L*sin(omega)+a*cos(lambda)/cos(omega)+a*sin(lambda)/cos(omega),a*sin(lambda)-a*cos(lambda)];
node.coordinates_mat(10,:)=[0,L*sin(omega)+a*sin(lambda)/cos(omega),-a*cos(lambda)];

bar.E_vec(19:22)=barE;
bar.A_vec(19:22)=barA;
bar.node_ij_mat(19,:)=[6,9];
bar.node_ij_mat(20,:)=[9,10];
bar.node_ij_mat(21,:)=[5,10];
bar.node_ij_mat(22,:)=[5,9];

node.coordinates_mat(11,:)=[-L*cos(omega),0,0];
node.coordinates_mat(12,:)=[-L*cos(omega)-a*cos(lambda)/sin(omega),0,a*sin(lambda)];
node.coordinates_mat(13,:)=[-L*cos(omega)-a*cos(lambda)/sin(omega)-a*sin(lambda)/sin(omega),0,a*sin(lambda)-a*cos(lambda)];
node.coordinates_mat(14,:)=[-L*cos(omega)-a*sin(lambda)/sin(omega),0,-a*cos(lambda)];

bar.E_vec(23:27)=barE;
bar.A_vec(23:27)=barA;
bar.node_ij_mat(23,:)=[11,12];
bar.node_ij_mat(24,:)=[12,13];
bar.node_ij_mat(25,:)=[13,14];
bar.node_ij_mat(26,:)=[11,14];
bar.node_ij_mat(27,:)=[11,13];

bar.E_vec(28:31)=barE;
bar.A_vec(28:31)=barA;
bar.node_ij_mat(28,:)=[5,11];
bar.node_ij_mat(29,:)=[6,12];
bar.node_ij_mat(30,:)=[9,13];
bar.node_ij_mat(31,:)=[10,14];

bar.E_vec(32:35)=barE;
bar.A_vec(32:35)=barA;
bar.node_ij_mat(32,:)=[5,12];
bar.node_ij_mat(33,:)=[6,13];
bar.node_ij_mat(34,:)=[9,14];
bar.node_ij_mat(35,:)=[10,11];

% Third Link
node.coordinates_mat(15,:)=[-L*cos(omega),0,0];
node.coordinates_mat(16,:)=[-L*cos(omega)-a*sin(lambda)/sin(omega),0,-a*cos(lambda)];

bar.E_vec(36:39)=barE;
bar.A_vec(36:39)=barA;
bar.node_ij_mat(36,:)=[12,15];
bar.node_ij_mat(37,:)=[13,16];
bar.node_ij_mat(38,:)=[15,16];
bar.node_ij_mat(39,:)=[12,16];

node.coordinates_mat(17,:)=[0,-L*sin(omega),0];
node.coordinates_mat(18,:)=[0,-L*sin(omega)-a*cos(lambda)/cos(omega),a*sin(lambda)];
node.coordinates_mat(19,:)=[0,-L*sin(omega)-a*cos(lambda)/cos(omega)-a*sin(lambda)/cos(omega),a*sin(lambda)-a*cos(lambda)];
node.coordinates_mat(20,:)=[0,-L*sin(omega)-a*sin(lambda)/cos(omega),-a*cos(lambda)];

bar.E_vec(40:43)=barE;
bar.A_vec(40:43)=barA;
bar.node_ij_mat(40,:)=[15,17];
bar.node_ij_mat(41,:)=[12,18];
bar.node_ij_mat(42,:)=[13,19];
bar.node_ij_mat(43,:)=[16,20];

bar.E_vec(44:47)=barE;
bar.A_vec(44:47)=barA;
bar.node_ij_mat(44,:)=[15,18];
bar.node_ij_mat(45,:)=[12,19];
bar.node_ij_mat(46,:)=[13,20];
bar.node_ij_mat(47,:)=[16,17];

bar.E_vec(48:52)=barE;
bar.A_vec(48:52)=barA;
bar.node_ij_mat(48,:)=[17,18];
bar.node_ij_mat(49,:)=[18,19];
bar.node_ij_mat(50,:)=[19,20];
bar.node_ij_mat(51,:)=[17,20];
bar.node_ij_mat(52,:)=[17,19];

% Forth Link
node.coordinates_mat(21,:)=[0,-L*sin(omega)-a*cos(lambda)/cos(omega)-a*sin(lambda)/cos(omega),a*sin(lambda)-a*cos(lambda)];
node.coordinates_mat(22,:)=[0,-L*sin(omega)-a*sin(lambda)/cos(omega),-a*cos(lambda)];

bar.E_vec(53:56)=barE;
bar.A_vec(53:56)=barA;
bar.node_ij_mat(53,:)=[18,21];
bar.node_ij_mat(54,:)=[21,22];
bar.node_ij_mat(55,:)=[17,22];
bar.node_ij_mat(56,:)=[17,21];

node.coordinates_mat(23,:)=[L*cos(omega),0,0];
node.coordinates_mat(24,:)=[L*cos(omega)+a*sin(lambda)/sin(omega),0,-a*cos(lambda)];

bar.E_vec(57:60)=barE;
bar.A_vec(57:60)=barA;
bar.node_ij_mat(57,:)=[2,23];
bar.node_ij_mat(58,:)=[3,24];
bar.node_ij_mat(59,:)=[23,24];
bar.node_ij_mat(60,:)=[3,23];

bar.E_vec(61:64)=barE;
bar.A_vec(61:64)=barA;
bar.node_ij_mat(61,:)=[17,23];
bar.node_ij_mat(62,:)=[2,18];
bar.node_ij_mat(63,:)=[3,21];
bar.node_ij_mat(64,:)=[22,24];

bar.E_vec(65:68)=barE;
bar.A_vec(65:68)=barA;
bar.node_ij_mat(65,:)=[2,17];
bar.node_ij_mat(66,:)=[2,21];
bar.node_ij_mat(67,:)=[3,22];
bar.node_ij_mat(68,:)=[22,23];

bar.Initialize(node);

% Add Rotational Hinge
rotSpr.node_ijkl_mat(1,:)=[2,5,6,12];
rotSpr.node_ijkl_mat(2,:)=[6,12,13,18];
rotSpr.node_ijkl_mat(3,:)=[12,17,18,2];
rotSpr.node_ijkl_mat(4,:)=[18,2,3,6];

rotSpr.rot_spr_K_vec(1,1)=stiffFactor*sprStiff;
rotSpr.rot_spr_K_vec(2,1)=sprStiff;
rotSpr.rot_spr_K_vec(3,1)=sprStiff;
rotSpr.rot_spr_K_vec(4,1)=sprStiff;

%% Initialize assembly
assembly=Assembly_Linkage();
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_4N=rotSpr;

assembly.InitializeAssembly()

%% Plot for investigation
plots=Plot_Linkage();
plots.displayRange=0.1;
plots.displayRangeRatio=1;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber()
plots.Plot_Shape_BarNumber()
plots.Plot_Shape_SprNumber()


plots.panelConnection{1}=[1,2,6,5];
plots.panelConnection{2}=[2,3,7,6];
plots.panelConnection{3}=[3,4,8,7];
plots.panelConnection{4}=[4,1,5,8];

plots.panelConnection{5}=[5 6 12 11];
plots.panelConnection{6}=[6 9 13 12];
plots.panelConnection{7}=[9 10 14 13];
plots.panelConnection{8}=[5 10 14 11];

plots.panelConnection{9}=[12 13 19 18];
plots.panelConnection{10}=[13 16 20 19];
plots.panelConnection{11}=[16 15 17 20];
plots.panelConnection{12}=[15 12 18 17];

plots.panelConnection{13}=[17 18 2 23];
plots.panelConnection{14}=[17 22 24 23];
plots.panelConnection{15}=[22 21 3 24];
plots.panelConnection{16}=[21 18 2 3];


%% Find alpha angle
% This is for comparing with the analytical solution
z1=node.coordinates_mat(6,:)-node.coordinates_mat(5,:);
x1=node.coordinates_mat(6,:)-node.coordinates_mat(2,:);
z2=node.coordinates_mat(13,:)-node.coordinates_mat(12,:);
x2=node.coordinates_mat(12,:)-node.coordinates_mat(6,:);

z1=z1/norm(z1);
x1=x1/norm(x1);
y1=cross(x1,z1);

z2=z2/norm(z2);
x2=x2/norm(x2);
y2=cross(x2,z2);

acos(dot(z1,z2))


%% Setup the loading controller
sf=Solver_NR_Folding_4N;
sf.assembly=assembly;
sf.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;
         4,1,1,1];

sf.targetRot=rotSpr.theta_current_vec;
sf.targetRot(1)=sf.targetRot(1)+0.7*pi;


sf.increStep=200;
sf.tol=10^-6;
sf.iterMax=50;

Uhis=sf.Solve();

plots.displayRange=0.15;
plots.fileName='Bennett_Linkage.gif';
plots.Plot_DeformedShape(squeeze(Uhis(150,:,:)))
plots.Plot_DeformedHis(Uhis(1:10:end,:,:))

