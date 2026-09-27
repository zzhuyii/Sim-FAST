clear all
close all
clc
tic

%% Define Geometry
% Size
length=0.05;

% initial folding status
alpha=(130/180)*pi; 

%% Define assembly
assembly=Assembly_Foldable_Unit;
cst=Vec_Elements_CST;
rot_spr_4N_D=Vec_Elements_RotSprings_4N_Directional_Smooth;
node=Elements_Nodes;

assembly.cst=cst;
assembly.rot_spr_4N_D=rot_spr_4N_D;
assembly.node=node;

%% Nodes Define
node.coordinates_mat=[node.coordinates_mat;
    length 0 0;
    0 length 0;
    0 0 0;
    length*cos(alpha) 0 length*sin(alpha);];

%% Define Plotting Functions
plots=Plot_Foldable_Unit;
plots.assembly=assembly;
plots.displayRange=[-0.05; 0.05; -0.05; 0.05; -0.05; 0.05];
plots.viewAngle1=20;
plots.viewAngle2=20;
plots.holdTime=0.04;

plots.Plot_Shape_Node_Number;

%% CST Define
tri_ijk=[];
tri_ijk=[tri_ijk;
    1 2 3;
    2 3 4;];

cst.node_ijk_mat=tri_ijk;
cstNum=size(cst.node_ijk_mat,1);

cst.t_vec=10*ones(cstNum,1);
cst.E_vec=4*10^9*ones(cstNum,1);
cst.v_vec=0.2*ones(cstNum,1);

plots.Plot_Shape_CST_Number;

%% Define Directional Rotational Spring
rot_spr_4N_D.node_ijkl_mat=[
            rot_spr_4N_D.node_ijkl_mat;
            1 2 3 4;];

rot_dir_num=size(rot_spr_4N_D.node_ijkl_mat);
rot_dir_num=rot_dir_num(1);

% Directional rotational springs for fold lines
rot_spr_4N_D.rot_spr_K_vec=0.0567*ones(rot_dir_num,1);
% Stiffness increase factor for directional springs when folded in the
% opposite direction
rot_spr_4N_D.mv_factor_vec=10*ones(rot_dir_num,1);

% The Mountain Valley Assignment for the directional spring
rot_spr_4N_D.mv_vec = ones(rot_dir_num,1);
rot_spr_4N_D.mv_vec(1) = 0;

plots.Plot_Shape_DirectSpr_Number();

assembly.Initialize_Assembly();


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;
nr.iterMax=40;
nr.tol=10^-6;
nr.increStep=1;

nr.supp=[
    1 1 1 1;    
    2 1 1 1;
    3 1 1 1;];

% This factor reduce the Newton step. 
% This is usually called damped NR solver.
nr.dampFactor=0.5;

% Total Steps
Uhis=[];

% Batch calculation settings
force_min = 0;
force_max = 5;
num_force = 400;

% Target nodal forces
nodal_force_vec = linspace(force_min, force_max, num_force);
nodeNum = size(node.coordinates_mat, 1);

% Preallocate result arrays
moment_vec    = zeros(num_force, 1);
theta_deg_vec = zeros(num_force, 1);

for i = 1:num_force
    nodal_force_magnitude = nodal_force_vec(i);
    

    % Load apply on folding
    nodeNum=size(node.coordinates_mat,1);
    nr.load = zeros(nodeNum, 4);
    nr.load(:,1) = (1:nodeNum)';
    nr.load(4, 2:4) = (nodal_force_magnitude)*[0,0,-1];

    Uhis(i,:,:)=squeeze(nr.Solve());
    xcurrent=node.coordinates_mat+node.current_U_mat;

    % Moment calculation
    x_distance = abs(xcurrent(4,1));
    moment_vec(i) = nodal_force_magnitude * x_distance;

    % Angle measurement
    p1 = xcurrent(1,:);
    p2 = xcurrent(2,:);
    p3 = xcurrent(3,:);
    p4 = xcurrent(4,:);

    a_m_n_1 = cross(p3 - p1, p2 - p1);
    a_m_n_2 = cross(p2 - p4, p3 - p4);
    cos_theta = dot(a_m_n_1, a_m_n_2) / (norm(a_m_n_1) * norm(a_m_n_2));
    cos_theta = max(-1,min(1,cos_theta));
    theta_rad = acos(cos_theta);
    theta_deg_vec(i) = rad2deg(theta_rad);


end
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)))
% plots.Plot_Deformed_Shape(U_final(:,:,end));

%% Plot Moment vs. Folding Angle
figure;

plot(theta_deg_vec, moment_vec, 'o-', ...
    'LineWidth', 2, ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', 'b');

xlabel('Folding Angle, \theta (degrees)');
ylabel('Moment (N·m)');
title('Moment vs. Folding Angle');

grid on;
box on;
