clear all
close all
clc
tic

%% Define Geometry

% Size of Skeleton
w=0.2;
h=0.2;
d=0.2;
t=0.02; 

% initial folding status
alpha=(89/180)*pi;  

% Size of Actuator
l_x=0.08; 
l_y=0.08;

% Layers amount
layer_num = 6; 

% each layer's angle
square_angle = (1/layer_num) * (pi-2*alpha); 

% Unit Number
N=6; 

% Target pressure to be applied
pressure=(144)*1000;

External_load_mass=0; % kg
% Target load on the structure
load_on_structure=(External_load_mass*9.81)/4;


%% Define assembly
assembly=Assembly_Foldable_Unit;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;
rot_spr_4N_D=Vec_Elements_RotSprings_4N_Directional_Smooth;
zlspr=Vec_Elements_Zero_L_Spring;
node=Elements_Nodes;

assembly.cst=cst;
assembly.node=node;
assembly.rot_spr_4N=rot_spr_4N;
assembly.rot_spr_4N_D=rot_spr_4N_D;
assembly.zlspr=zlspr;


%% Nodes Define

% unit 1
i=1;
node.coordinates_mat=[node.coordinates_mat;
    -((i-1)*(w*cos(alpha)+2*t)+0), d, 0;
    -((i-1)*(w*cos(alpha)+2*t)+0), 0, 0;
    -((i-1)*(w*cos(alpha)+2*t)+0), d, h;
    -((i-1)*(w*cos(alpha)+2*t)+0), 0, h;
    -((i-1)*(w*cos(alpha)+2*t)+t), d, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t), 0, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t), d, h;
    -((i-1)*(w*cos(alpha)+2*t)+t), 0, h; % 8 nodes

    -((i-1)*(w*cos(alpha)+2*t)+t), d, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t), d/2, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t), 0, 0; % 11 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d, (w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d-((d-l_y)/2), (w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d-((d-l_y)/2)-l_y, (w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), 0, (w/2-l_x)*sin(alpha); % 15 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d, (w/2)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d-((d-l_y)/2), (w/2)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d/2, (w/2)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), (d-l_y)/2, (w/2)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), 0, (w/2)*sin(alpha); % 20 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d, (w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d-((d-l_y)/2), (w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d-((d-l_y)/2)-l_y, (w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), 0, (w/2-l_x)*sin(alpha); % 24 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), d, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), d/2, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), 0, 0; % 27 nodes

    -((i-1)*(w*cos(alpha)+2*t)+t), d, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t), d, h/2;
    -((i-1)*(w*cos(alpha)+2*t)+t), d, h; % 30 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d+(w/2-l_x)*sin(alpha), 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d+(w/2-l_x)*sin(alpha), h-((h-l_y)/2)-l_y;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d+(w/2-l_x)*sin(alpha), h-((h-l_y)/2);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d+(w/2-l_x)*sin(alpha), h; % 34 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d+(w/2)*sin(alpha), 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d+(w/2)*sin(alpha), h-((h-l_y)/2)-l_y;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d+(w/2)*sin(alpha), h/2;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d+(w/2)*sin(alpha), h-((h-l_y)/2);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d+(w/2)*sin(alpha), h; % 39 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d+(w/2-l_x)*sin(alpha), 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d+(w/2-l_x)*sin(alpha), h-((h-l_y)/2)-l_y;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d+(w/2-l_x)*sin(alpha), h-((h-l_y)/2);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d+(w/2-l_x)*sin(alpha), h; % 43 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), d, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), d, h/2;
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), d, h; % 46 nodes

    -((i-1)*(w*cos(alpha)+2*t)+t), d, h;
    -((i-1)*(w*cos(alpha)+2*t)+t), d/2, h;
    -((i-1)*(w*cos(alpha)+2*t)+t), 0, h; % 49 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d, h-(w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d-((d-l_y)/2), h-(w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), d-((d-l_y)/2)-l_y, h-(w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), 0, h-(w/2-l_x)*sin(alpha); % 53 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d, h-(w/2)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d-((d-l_y)/2), h-(w/2)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), d/2, h-(w/2)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), (d-l_y)/2, h-(w/2)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), 0, h-(w/2)*sin(alpha); % 58 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d, h-(w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d-((d-l_y)/2), h-(w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), d-((d-l_y)/2)-l_y, h-(w/2-l_x)*sin(alpha);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), 0, h-(w/2-l_x)*sin(alpha); % 62 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), d, h;
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), d/2, h;
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), 0, h; % 65 nodes

    -((i-1)*(w*cos(alpha)+2*t)+t), 0, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t), 0, h/2;
    -((i-1)*(w*cos(alpha)+2*t)+t), 0, h; % 68 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), 0-(w/2-l_x)*sin(alpha), 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), 0-(w/2-l_x)*sin(alpha), h-((h-l_y)/2)-l_y;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), 0-(w/2-l_x)*sin(alpha), h-((h-l_y)/2);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2-l_x)*cos(alpha)), 0-(w/2-l_x)*sin(alpha), h; % 72 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), 0-(w/2)*sin(alpha), 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), 0-(w/2)*sin(alpha), h-((h-l_y)/2)-l_y;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), 0-(w/2)*sin(alpha), h/2;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), 0-(w/2)*sin(alpha), h-((h-l_y)/2);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2)*cos(alpha)), 0-(w/2)*sin(alpha), h; % 77 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), 0-(w/2-l_x)*sin(alpha), 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), 0-(w/2-l_x)*sin(alpha), h-((h-l_y)/2)-l_y;
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), 0-(w/2-l_x)*sin(alpha), h-((h-l_y)/2);
    -((i-1)*(w*cos(alpha)+2*t)+t+(w/2+l_x)*cos(alpha)), 0-(w/2-l_x)*sin(alpha), h; % 81 nodes
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), 0, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), 0, h/2;
    -((i-1)*(w*cos(alpha)+2*t)+t+w*cos(alpha)), 0, h; % 84 nodes

    -((i-1)*(w*cos(alpha)+2*t)+0+t+w*cos(alpha)), d, 0;
    -((i-1)*(w*cos(alpha)+2*t)+0+t+w*cos(alpha)), 0, 0;
    -((i-1)*(w*cos(alpha)+2*t)+0+t+w*cos(alpha)), d, h;
    -((i-1)*(w*cos(alpha)+2*t)+0+t+w*cos(alpha)), 0, h;
    -((i-1)*(w*cos(alpha)+2*t)+t+t+w*cos(alpha)), d, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+t+w*cos(alpha)), 0, 0;
    -((i-1)*(w*cos(alpha)+2*t)+t+t+w*cos(alpha)), d, h;
    -((i-1)*(w*cos(alpha)+2*t)+t+t+w*cos(alpha)), 0, h;]; % 92 nodes
% Rotate all 92 nodes around the node 1–2 line
rotation_angle = deg2rad(75);  % Use -75 for the opposite direction

R_y = [
     cos(rotation_angle), 0, sin(rotation_angle);
     0,                   1, 0;
    -sin(rotation_angle), 0, cos(rotation_angle)
];

% Current unit: the most recently added 92 nodes
node_range = (size(node.coordinates_mat,1)-91):size(node.coordinates_mat,1);

% Node 1 is a point on the node 1–2 rotation axis
axis_point = node.coordinates_mat(node_range(1),:);

% Rotate all 92 nodes
coordinates = node.coordinates_mat(node_range,:);
coordinates = coordinates - axis_point;
coordinates = (R_y * coordinates')';
coordinates = coordinates + axis_point;

node.coordinates_mat(node_range,:) = coordinates;

% Copy the rotated unit around a new y-direction axis

nodes_per_unit = 92;

% Save the current 92-node unit after its 75-degree rotation
original_unit = node.coordinates_mat(node_range,:);

% Original node 1-2 axis
original_axis_point = original_unit(1,:);

% Distance from the node 1-2 line to the new rotation axis
rotation_radius = ((norm(node.coordinates_mat(1,:)-node.coordinates_mat(89,:)))/2) / sind(15);

% New axis is parallel to the y-axis and shifted in the negative x-direction
new_axis_point = original_axis_point + [-rotation_radius, 0, 0];

% Number of additional copies
num_additional_units = 5;   % Change this value as needed

% Angular spacing between adjacent units
angle_increment = deg2rad(-30);

for copy_id = 1:num_additional_units

    rotation_angle_copy = copy_id * angle_increment;

    R_y_copy = [
         cos(rotation_angle_copy), 0, sin(rotation_angle_copy);
         0,                        1, 0;
        -sin(rotation_angle_copy), 0, cos(rotation_angle_copy)
    ];

    % Start every copy from the original 75-degree-rotated unit
    copied_unit = original_unit;

    % Rotate about the new y-direction axis
    copied_unit = copied_unit - new_axis_point;
    copied_unit = (R_y_copy * copied_unit')';
    copied_unit = copied_unit + new_axis_point;

    % Add the copied unit without deleting the original unit
    node.coordinates_mat = [
        node.coordinates_mat;
        copied_unit
    ];
end


% Connection nodes
node.coordinates_mat = [
        node.coordinates_mat;
        node.coordinates_mat(1,:); % 553
        node.coordinates_mat(2,:);
        node.coordinates_mat(3,:);
        node.coordinates_mat(4,:);
        node.coordinates_mat(3,1),node.coordinates_mat(3,2),0;
        node.coordinates_mat(4,1),node.coordinates_mat(4,2),0;];

node.coordinates_mat = [
        node.coordinates_mat;
        node.coordinates_mat(89,:); % 559
        node.coordinates_mat(90,:);
        node.coordinates_mat(91,:);
        node.coordinates_mat(92,:);
        node.coordinates_mat(95,:);
        node.coordinates_mat(96,:);];

node.coordinates_mat = [
        node.coordinates_mat;
        node.coordinates_mat(181,:); % 565
        node.coordinates_mat(182,:);
        node.coordinates_mat(183,:);
        node.coordinates_mat(184,:);
        node.coordinates_mat(187,:);
        node.coordinates_mat(188,:);];

node.coordinates_mat = [
        node.coordinates_mat;
        node.coordinates_mat(273,:); % 571
        node.coordinates_mat(274,:);
        node.coordinates_mat(275,:);
        node.coordinates_mat(276,:);
        node.coordinates_mat(279,:);
        node.coordinates_mat(280,:);];

node.coordinates_mat = [
        node.coordinates_mat;
        node.coordinates_mat(365,:); % 577
        node.coordinates_mat(366,:);
        node.coordinates_mat(367,:);
        node.coordinates_mat(368,:);
        node.coordinates_mat(371,:);
        node.coordinates_mat(372,:);];

node.coordinates_mat = [
        node.coordinates_mat;
        node.coordinates_mat(457,:); % 583
        node.coordinates_mat(458,:);
        node.coordinates_mat(459,:);
        node.coordinates_mat(460,:);
        node.coordinates_mat(463,:);
        node.coordinates_mat(464,:);];

node.coordinates_mat = [
        node.coordinates_mat;
        node.coordinates_mat(549,:); % 589
        node.coordinates_mat(550,:);
        node.coordinates_mat(551,:);
        node.coordinates_mat(552,:);
        node.coordinates_mat(551,1),node.coordinates_mat(551,2),0;
        node.coordinates_mat(552,1),node.coordinates_mat(552,2),0;];

%% Define Plotting Functions
plots=Plot_Foldable_Unit;
plots.assembly=assembly;
plots.displayRange=[-1.2; 0.4; -0.5; 0.5; -0.1; 1];
plots.viewAngle1=20;
plots.viewAngle2=20;
plots.holdTime=0.04;

plots.Plot_Shape_Node_Number;


%% CST Define
tri_ijk=[];
tri_direction=[];
% n_e_u = size(node.coordinates_mat, 1)/N; % each unit nodes amount
n_e_u = 92;

for i=1:N
    % Skeleton
     tri_ijk=[tri_ijk;
        n_e_u*(i-1)+1    n_e_u*(i-1)+5    n_e_u*(i-1)+6;
        n_e_u*(i-1)+1    n_e_u*(i-1)+2    n_e_u*(i-1)+6;
        n_e_u*(i-1)+1    n_e_u*(i-1)+5    n_e_u*(i-1)+7;
        n_e_u*(i-1)+1    n_e_u*(i-1)+3    n_e_u*(i-1)+7;
        n_e_u*(i-1)+3    n_e_u*(i-1)+7    n_e_u*(i-1)+8;
        n_e_u*(i-1)+3    n_e_u*(i-1)+4    n_e_u*(i-1)+8;
        n_e_u*(i-1)+4    n_e_u*(i-1)+8    n_e_u*(i-1)+6;
        n_e_u*(i-1)+4    n_e_u*(i-1)+2    n_e_u*(i-1)+6;];  % 8

    for j=1:4
         tri_ijk=[tri_ijk;
             n_e_u*(i-1)+(j-1)*19+9    n_e_u*(i-1)+(j-1)*19+10    n_e_u*(i-1)+(j-1)*19+13;
             n_e_u*(i-1)+(j-1)*19+10    n_e_u*(i-1)+(j-1)*19+11    n_e_u*(i-1)+(j-1)*19+14;
             n_e_u*(i-1)+(j-1)*19+9    n_e_u*(i-1)+(j-1)*19+12    n_e_u*(i-1)+(j-1)*19+13;
             n_e_u*(i-1)+(j-1)*19+10    n_e_u*(i-1)+(j-1)*19+13    n_e_u*(i-1)+(j-1)*19+14;
             n_e_u*(i-1)+(j-1)*19+11    n_e_u*(i-1)+(j-1)*19+14    n_e_u*(i-1)+(j-1)*19+15;
             n_e_u*(i-1)+(j-1)*19+12    n_e_u*(i-1)+(j-1)*19+13    n_e_u*(i-1)+(j-1)*19+16;
             n_e_u*(i-1)+(j-1)*19+13    n_e_u*(i-1)+(j-1)*19+16    n_e_u*(i-1)+(j-1)*19+17;
             n_e_u*(i-1)+(j-1)*19+13    n_e_u*(i-1)+(j-1)*19+17    n_e_u*(i-1)+(j-1)*19+18;
             n_e_u*(i-1)+(j-1)*19+13    n_e_u*(i-1)+(j-1)*19+14    n_e_u*(i-1)+(j-1)*19+18;
             n_e_u*(i-1)+(j-1)*19+14    n_e_u*(i-1)+(j-1)*19+18    n_e_u*(i-1)+(j-1)*19+19;
             n_e_u*(i-1)+(j-1)*19+14    n_e_u*(i-1)+(j-1)*19+19    n_e_u*(i-1)+(j-1)*19+20;
             n_e_u*(i-1)+(j-1)*19+14    n_e_u*(i-1)+(j-1)*19+15    n_e_u*(i-1)+(j-1)*19+20;

             n_e_u*(i-1)+(j-1)*19+16    n_e_u*(i-1)+(j-1)*19+21    n_e_u*(i-1)+(j-1)*19+22;
             n_e_u*(i-1)+(j-1)*19+16    n_e_u*(i-1)+(j-1)*19+17    n_e_u*(i-1)+(j-1)*19+22;
             n_e_u*(i-1)+(j-1)*19+17    n_e_u*(i-1)+(j-1)*19+18    n_e_u*(i-1)+(j-1)*19+22;
             n_e_u*(i-1)+(j-1)*19+18    n_e_u*(i-1)+(j-1)*19+22    n_e_u*(i-1)+(j-1)*19+23;
             n_e_u*(i-1)+(j-1)*19+18    n_e_u*(i-1)+(j-1)*19+19    n_e_u*(i-1)+(j-1)*19+23;
             n_e_u*(i-1)+(j-1)*19+19    n_e_u*(i-1)+(j-1)*19+20    n_e_u*(i-1)+(j-1)*19+23;
             n_e_u*(i-1)+(j-1)*19+20    n_e_u*(i-1)+(j-1)*19+23    n_e_u*(i-1)+(j-1)*19+24;
             n_e_u*(i-1)+(j-1)*19+21    n_e_u*(i-1)+(j-1)*19+22    n_e_u*(i-1)+(j-1)*19+25;
             n_e_u*(i-1)+(j-1)*19+22    n_e_u*(i-1)+(j-1)*19+25    n_e_u*(i-1)+(j-1)*19+26;
             n_e_u*(i-1)+(j-1)*19+22    n_e_u*(i-1)+(j-1)*19+23    n_e_u*(i-1)+(j-1)*19+26;
             n_e_u*(i-1)+(j-1)*19+23    n_e_u*(i-1)+(j-1)*19+26    n_e_u*(i-1)+(j-1)*19+27;
             n_e_u*(i-1)+(j-1)*19+23    n_e_u*(i-1)+(j-1)*19+24    n_e_u*(i-1)+(j-1)*19+27;];    % 104
    end
    tri_ijk=[tri_ijk;
        n_e_u*(i-1)+85    n_e_u*(i-1)+89    n_e_u*(i-1)+90;
        n_e_u*(i-1)+85    n_e_u*(i-1)+86    n_e_u*(i-1)+90;
        n_e_u*(i-1)+85    n_e_u*(i-1)+89    n_e_u*(i-1)+91;
        n_e_u*(i-1)+85    n_e_u*(i-1)+87    n_e_u*(i-1)+91;
        n_e_u*(i-1)+87    n_e_u*(i-1)+91    n_e_u*(i-1)+92;
        n_e_u*(i-1)+87    n_e_u*(i-1)+88    n_e_u*(i-1)+92;
        n_e_u*(i-1)+88    n_e_u*(i-1)+92    n_e_u*(i-1)+90;
        n_e_u*(i-1)+88    n_e_u*(i-1)+86    n_e_u*(i-1)+90;];    % 112
end

for r=1:7
tri_ijk=[tri_ijk;
    553+(r-1)*6 555+(r-1)*6 557+(r-1)*6;
    554+(r-1)*6 556+(r-1)*6 558+(r-1)*6;
    553+(r-1)*6 555+(r-1)*6 556+(r-1)*6;
    553+(r-1)*6 554+(r-1)*6 556+(r-1)*6;
    553+(r-1)*6 557+(r-1)*6 558+(r-1)*6;
    553+(r-1)*6 554+(r-1)*6 558+(r-1)*6;
    557+(r-1)*6 555+(r-1)*6 556+(r-1)*6;
    557+(r-1)*6 558+(r-1)*6 556+(r-1)*6;];
end

cst.node_ijk_mat=tri_ijk;
cstNum=size(cst.node_ijk_mat,1);

cst.t_vec=0.006*ones(cstNum,1);
cst.E_vec=4*10^9*ones(cstNum,1);
cst.v_vec=0.2*ones(cstNum,1);

plots.Plot_Shape_CST_Number;



%% Define Rotational Spring
for i=1:N
    % Skeleton
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        n_e_u*(i-1)+6    n_e_u*(i-1)+1    n_e_u*(i-1)+5    n_e_u*(i-1)+7;
        n_e_u*(i-1)+1    n_e_u*(i-1)+3    n_e_u*(i-1)+7    n_e_u*(i-1)+8;
        n_e_u*(i-1)+3    n_e_u*(i-1)+4    n_e_u*(i-1)+8    n_e_u*(i-1)+6;
        n_e_u*(i-1)+4    n_e_u*(i-1)+2    n_e_u*(i-1)+6    n_e_u*(i-1)+1;
        n_e_u*(i-1)+5    n_e_u*(i-1)+1    n_e_u*(i-1)+7    n_e_u*(i-1)+3;
        n_e_u*(i-1)+7    n_e_u*(i-1)+3    n_e_u*(i-1)+8    n_e_u*(i-1)+4;
        n_e_u*(i-1)+8    n_e_u*(i-1)+4    n_e_u*(i-1)+6    n_e_u*(i-1)+2;
        n_e_u*(i-1)+5    n_e_u*(i-1)+1    n_e_u*(i-1)+6    n_e_u*(i-1)+2;];

        for j=1:4
        rot_spr_4N.node_ijkl_mat=[
            rot_spr_4N.node_ijkl_mat;
            n_e_u*(i-1)+19*(j-1)+9    n_e_u*(i-1)+19*(j-1)+10    n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+14;
            n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+10    n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+11;
            n_e_u*(i-1)+19*(j-1)+10    n_e_u*(i-1)+19*(j-1)+11    n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+15;
            n_e_u*(i-1)+19*(j-1)+11    n_e_u*(i-1)+19*(j-1)+15    n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+20;
            n_e_u*(i-1)+19*(j-1)+15    n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+20    n_e_u*(i-1)+19*(j-1)+19;
            n_e_u*(i-1)+19*(j-1)+20    n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+19    n_e_u*(i-1)+19*(j-1)+18;
            n_e_u*(i-1)+19*(j-1)+19    n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+13;
            n_e_u*(i-1)+19*(j-1)+10    n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+18;
            n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+17;
            n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+17    n_e_u*(i-1)+19*(j-1)+16;
            n_e_u*(i-1)+19*(j-1)+17    n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+16    n_e_u*(i-1)+19*(j-1)+12;
            n_e_u*(i-1)+19*(j-1)+16    n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+12    n_e_u*(i-1)+19*(j-1)+9;
            n_e_u*(i-1)+19*(j-1)+10    n_e_u*(i-1)+19*(j-1)+9    n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+12;

            n_e_u*(i-1)+19*(j-1)+16    n_e_u*(i-1)+19*(j-1)+17    n_e_u*(i-1)+19*(j-1)+22    n_e_u*(i-1)+19*(j-1)+18;
            n_e_u*(i-1)+19*(j-1)+17    n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+22    n_e_u*(i-1)+19*(j-1)+23;
            n_e_u*(i-1)+19*(j-1)+22    n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+23    n_e_u*(i-1)+19*(j-1)+19;
            n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+19    n_e_u*(i-1)+19*(j-1)+23    n_e_u*(i-1)+19*(j-1)+20;
            n_e_u*(i-1)+19*(j-1)+19    n_e_u*(i-1)+19*(j-1)+20    n_e_u*(i-1)+19*(j-1)+23    n_e_u*(i-1)+19*(j-1)+24;
            n_e_u*(i-1)+19*(j-1)+20    n_e_u*(i-1)+19*(j-1)+24    n_e_u*(i-1)+19*(j-1)+23    n_e_u*(i-1)+19*(j-1)+27;
            n_e_u*(i-1)+19*(j-1)+24    n_e_u*(i-1)+19*(j-1)+23    n_e_u*(i-1)+19*(j-1)+27    n_e_u*(i-1)+19*(j-1)+26;
            n_e_u*(i-1)+19*(j-1)+27    n_e_u*(i-1)+19*(j-1)+23    n_e_u*(i-1)+19*(j-1)+26    n_e_u*(i-1)+19*(j-1)+22;
            n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+22    n_e_u*(i-1)+19*(j-1)+23    n_e_u*(i-1)+19*(j-1)+26;
            n_e_u*(i-1)+19*(j-1)+23    n_e_u*(i-1)+19*(j-1)+22    n_e_u*(i-1)+19*(j-1)+26    n_e_u*(i-1)+19*(j-1)+25;
            n_e_u*(i-1)+19*(j-1)+26    n_e_u*(i-1)+19*(j-1)+22    n_e_u*(i-1)+19*(j-1)+25    n_e_u*(i-1)+19*(j-1)+21;
            n_e_u*(i-1)+19*(j-1)+25    n_e_u*(i-1)+19*(j-1)+22    n_e_u*(i-1)+19*(j-1)+21    n_e_u*(i-1)+19*(j-1)+16;
            n_e_u*(i-1)+19*(j-1)+21    n_e_u*(i-1)+19*(j-1)+16    n_e_u*(i-1)+19*(j-1)+22    n_e_u*(i-1)+19*(j-1)+17;];
        end
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        n_e_u*(i-1)+90    n_e_u*(i-1)+85    n_e_u*(i-1)+89    n_e_u*(i-1)+91;
        n_e_u*(i-1)+85    n_e_u*(i-1)+87    n_e_u*(i-1)+91    n_e_u*(i-1)+92;
        n_e_u*(i-1)+87    n_e_u*(i-1)+88    n_e_u*(i-1)+92    n_e_u*(i-1)+90;
        n_e_u*(i-1)+88    n_e_u*(i-1)+86    n_e_u*(i-1)+90    n_e_u*(i-1)+85;
        n_e_u*(i-1)+89    n_e_u*(i-1)+85    n_e_u*(i-1)+91    n_e_u*(i-1)+87;
        n_e_u*(i-1)+91    n_e_u*(i-1)+87    n_e_u*(i-1)+92    n_e_u*(i-1)+88;
        n_e_u*(i-1)+92    n_e_u*(i-1)+88    n_e_u*(i-1)+90    n_e_u*(i-1)+86;
        n_e_u*(i-1)+86    n_e_u*(i-1)+85    n_e_u*(i-1)+90    n_e_u*(i-1)+89;]; % 120

        for j=1:4
        rot_spr_4N_D.node_ijkl_mat=[
            rot_spr_4N_D.node_ijkl_mat;
            n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+16    n_e_u*(i-1)+19*(j-1)+17    n_e_u*(i-1)+19*(j-1)+22;
            n_e_u*(i-1)+19*(j-1)+13    n_e_u*(i-1)+19*(j-1)+17    n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+22;
            n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+18    n_e_u*(i-1)+19*(j-1)+19    n_e_u*(i-1)+19*(j-1)+23;
            n_e_u*(i-1)+19*(j-1)+14    n_e_u*(i-1)+19*(j-1)+19    n_e_u*(i-1)+19*(j-1)+20    n_e_u*(i-1)+19*(j-1)+23;]; % 136     Folding line 
        end
end

for r=1:7
rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        555+(r-1)*6 553+(r-1)*6 557+(r-1)*6 558+(r-1)*6;
        557+(r-1)*6 553+(r-1)*6 558+(r-1)*6 554+(r-1)*6;
        553+(r-1)*6 554+(r-1)*6 558+(r-1)*6 556+(r-1)*6;
        558+(r-1)*6 556+(r-1)*6 554+(r-1)*6 553+(r-1)*6;
        554+(r-1)*6 553+(r-1)*6 556+(r-1)*6 555+(r-1)*6;
        556+(r-1)*6 553+(r-1)*6 555+(r-1)*6 557+(r-1)*6;
        556+(r-1)*6 554+(r-1)*6 553+(r-1)*6 558+(r-1)*6;
        553+(r-1)*6 555+(r-1)*6 556+(r-1)*6 557+(r-1)*6;
        555+(r-1)*6 556+(r-1)*6 557+(r-1)*6 558+(r-1)*6;
        556+(r-1)*6 557+(r-1)*6 558+(r-1)*6 553+(r-1)*6;
        553+(r-1)*6 555+(r-1)*6 557+(r-1)*6 556+(r-1)*6;
        557+(r-1)*6 556+(r-1)*6 558+(r-1)*6 554+(r-1)*6;];
end

rotNum=size(rot_spr_4N.node_ijkl_mat);
rotNum=rotNum(1);
rot_dir_num=size(rot_spr_4N_D.node_ijkl_mat);
rot_dir_num=rot_dir_num(1);

% Stiff rotational springs for panels
rot_spr_4N.rot_spr_K_vec=5000*ones(rotNum,1); 

% Directional rotational springs for fold lines
rot_spr_4N_D.rot_spr_K_vec=0.015*ones(rot_dir_num,1); 

% Stiffness increase factor for directional springs when folded in the
% opposite direction
rot_spr_4N_D.mv_factor_vec = 500 * ones(rot_dir_num, 1);

% % Springs 33–64 use a larger increase factor
% rot_spr_4N_D.mv_factor_vec(33:64) = 1000;

% The Mountain Valley Assignment for the directional spring
rot_spr_4N_D.mv_vec = ones(rot_dir_num,1);

for i = 1:N
    idx = (i-1)*16 + [13:16, 9:12];   
    rot_spr_4N_D.mv_vec(idx) = 0;
end

plots.Plot_Shape_Spr_Number();
plots.Plot_Shape_DirectSpr_Number();


%% Define the connectors
zlsprStiff=10000000;
for i=1:N 
    zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+5 (i-1)*n_e_u+9; (i-1)*n_e_u+5 (i-1)*n_e_u+28];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+7 (i-1)*n_e_u+30; (i-1)*n_e_u+7 (i-1)*n_e_u+47];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+8 (i-1)*n_e_u+49; (i-1)*n_e_u+8 (i-1)*n_e_u+68];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+6 (i-1)*n_e_u+66; (i-1)*n_e_u+6 (i-1)*n_e_u+11];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+85 (i-1)*n_e_u+25; (i-1)*n_e_u+85 (i-1)*n_e_u+44];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+87 (i-1)*n_e_u+46; (i-1)*n_e_u+87 (i-1)*n_e_u+63];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+88 (i-1)*n_e_u+65; (i-1)*n_e_u+88 (i-1)*n_e_u+84];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+86 (i-1)*n_e_u+82; (i-1)*n_e_u+86 (i-1)*n_e_u+27];
end

zlspr.node_ij_mat=[zlspr.node_ij_mat;553 1];
zlspr.node_ij_mat=[zlspr.node_ij_mat;554 2];
zlspr.node_ij_mat=[zlspr.node_ij_mat;555 3];
zlspr.node_ij_mat=[zlspr.node_ij_mat;556 4];

for r=1:5
zlspr.node_ij_mat=[zlspr.node_ij_mat;559+(r-1)*6 89+(r-1)*92];
zlspr.node_ij_mat=[zlspr.node_ij_mat;560+(r-1)*6 90+(r-1)*92];
zlspr.node_ij_mat=[zlspr.node_ij_mat;561+(r-1)*6 91+(r-1)*92];
zlspr.node_ij_mat=[zlspr.node_ij_mat;562+(r-1)*6 92+(r-1)*92];
zlspr.node_ij_mat=[zlspr.node_ij_mat;559+(r-1)*6 93+(r-1)*92];
zlspr.node_ij_mat=[zlspr.node_ij_mat;560+(r-1)*6 94+(r-1)*92];
zlspr.node_ij_mat=[zlspr.node_ij_mat;563+(r-1)*6 95+(r-1)*92];
zlspr.node_ij_mat=[zlspr.node_ij_mat;564+(r-1)*6 96+(r-1)*92];
end

zlspr.node_ij_mat=[zlspr.node_ij_mat;589 549];
zlspr.node_ij_mat=[zlspr.node_ij_mat;590 550];
zlspr.node_ij_mat=[zlspr.node_ij_mat;591 551];
zlspr.node_ij_mat=[zlspr.node_ij_mat;592 552];

zlsprNum=size(zlspr.node_ij_mat,1);
zlspr.k_vec=zlsprStiff*ones(zlsprNum,1);

plots.Plot_Shape_ZLsprNumber();

assembly.Initialize_Assembly();

%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;
nr.iterMax=15;
nr.tol=10^-3;
nr.increStep=1;

nr.supp=[
    553 1 1 1;    
    554 0 1 1;
    557 0 1 1;
    558 0 1 1;
    589 0 1 1;
    590 0 1 1;
    593 0 1 1;
    594 0 1 1;];

% This factor reduce the Newton step. 
% This is usually called damped NR solver.
nr.dampFactor=0.5;

% Total Steps
step=200;

% Equivalent nodal force magnitude
equal_nodal_force_magnitude=1/4*(pressure*l_x*l_y);


Uhis=[];
for k=1:step

    xcurrent=node.coordinates_mat+node.current_U_mat;
    
    % Actuator pressure direction
    equal_nodal_force_direct = zeros(8*N,4);
    equal_nodal_force_direct(:,1) = (1:8*N)';
    for i=1:N
        
        n1=cross((xcurrent(19+n_e_u*(i-1),:)-xcurrent(17+n_e_u*(i-1),:)),(xcurrent(13+n_e_u*(i-1),:)-xcurrent(17+n_e_u*(i-1),:)));
        n2=cross((xcurrent(22+n_e_u*(i-1),:)-xcurrent(17+n_e_u*(i-1),:)),(xcurrent(19+n_e_u*(i-1),:)-xcurrent(17+n_e_u*(i-1),:)));

        n3=cross((xcurrent(38+n_e_u*(i-1),:)-xcurrent(36+n_e_u*(i-1),:)),(xcurrent(32+n_e_u*(i-1),:)-xcurrent(36+n_e_u*(i-1),:)));
        n4=cross((xcurrent(41+n_e_u*(i-1),:)-xcurrent(36+n_e_u*(i-1),:)),(xcurrent(38+n_e_u*(i-1),:)-xcurrent(36+n_e_u*(i-1),:)));

        n5=cross((xcurrent(51+n_e_u*(i-1),:)-xcurrent(55+n_e_u*(i-1),:)),(xcurrent(57+n_e_u*(i-1),:)-xcurrent(55+n_e_u*(i-1),:)));
        n6=cross((xcurrent(57+n_e_u*(i-1),:)-xcurrent(55+n_e_u*(i-1),:)),(xcurrent(60+n_e_u*(i-1),:)-xcurrent(55+n_e_u*(i-1),:)));

        n7=cross((xcurrent(70+n_e_u*(i-1),:)-xcurrent(74+n_e_u*(i-1),:)),(xcurrent(76+n_e_u*(i-1),:)-xcurrent(74+n_e_u*(i-1),:)));
        n8=cross((xcurrent(76+n_e_u*(i-1),:)-xcurrent(74+n_e_u*(i-1),:)),(xcurrent(79+n_e_u*(i-1),:)-xcurrent(74+n_e_u*(i-1),:)));

        n1=safeUnitNormal(n1,'n1',k);
        n2=safeUnitNormal(n2,'n2',k);
        n3=safeUnitNormal(n3,'n3',k);
        n4=safeUnitNormal(n4,'n4',k);
        n5=safeUnitNormal(n5,'n5',k);
        n6=safeUnitNormal(n6,'n6',k);
        n7=safeUnitNormal(n7,'n7',k);
        n8=safeUnitNormal(n8,'n8',k);

        equal_nodal_force_direct(1+8*(i-1),2:4)=n1;
        equal_nodal_force_direct(2+8*(i-1),2:4)=n2;
        equal_nodal_force_direct(3+8*(i-1),2:4)=n3;
        equal_nodal_force_direct(4+8*(i-1),2:4)=n4;
        equal_nodal_force_direct(5+8*(i-1),2:4)=n5;
        equal_nodal_force_direct(6+8*(i-1),2:4)=n6;
        equal_nodal_force_direct(7+8*(i-1),2:4)=n7;
        equal_nodal_force_direct(8+8*(i-1),2:4)=n8;
    end

    % Nodal force apply
    nodeNum=size(node.coordinates_mat,1);
    nr.load = zeros(nodeNum, 4);
    nr.load(:,1) = (1:nodeNum)';
    for i=1:N
        for j=1:4
        nr.load(13+19*(j-1)+n_e_u*(i-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*equal_nodal_force_direct(1+2*(j-1)+8*(i-1),2:4);
        nr.load(14+19*(j-1)+n_e_u*(i-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*equal_nodal_force_direct(1+2*(j-1)+8*(i-1),2:4);
        nr.load(17+19*(j-1)+n_e_u*(i-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*-equal_nodal_force_direct(1+2*(j-1)+8*(i-1),2:4)+((k)*(equal_nodal_force_magnitude/step))*-equal_nodal_force_direct(2+2*(j-1)+8*(i-1),2:4);
        nr.load(19+19*(j-1)+n_e_u*(i-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*-equal_nodal_force_direct(1+2*(j-1)+8*(i-1),2:4)+((k)*(equal_nodal_force_magnitude/step))*-equal_nodal_force_direct(2+2*(j-1)+8*(i-1),2:4);
        nr.load(22+19*(j-1)+n_e_u*(i-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*equal_nodal_force_direct(2+2*(j-1)+8*(i-1),2:4);
        nr.load(23+19*(j-1)+n_e_u*(i-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*equal_nodal_force_direct(2+2*(j-1)+8*(i-1),2:4);
        end
    end
    
    % % load put on the structure
    % nr.load(89, 2:4) = nr.load(89, 2:4) + [(k)*(load_on_structure/step) 0 0];
    % nr.load(90, 2:4) = nr.load(90, 2:4) + [(k)*(load_on_structure/step) 0 0];
    % nr.load(91, 2:4) = nr.load(91, 2:4) + [(k)*(load_on_structure/step) 0 0];
    % nr.load(92, 2:4) = nr.load(92, 2:4) + [(k)*(load_on_structure/step) 0 0];    

    Uhis(k,:,:)=squeeze(nr.Solve());

end

plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)))

plots.fileName='Demo_Arch.gif';
plots.Plot_Deformed_His(Uhis(1:1:end,:,:))


%% Angle measurement
% p31 = xcurrent(31,:);
% p35 = xcurrent(35,:);
% p36 = xcurrent(36,:);
% p40 = xcurrent(40,:);
% a_m_n_1 = cross(p35 - p31, p36 - p31);
% a_m_n_2 = cross(p36 - p35, p40 - p35);
% cos_theta = dot(a_m_n_1, a_m_n_2) / (norm(a_m_n_1) * norm(a_m_n_2));
% cos_theta = max(-1,min(1,cos_theta));
% theta_rad = acos(cos_theta);
% theta_deg = rad2deg(theta_rad); 


%% Force computation 
function n = safeUnitNormal(n,label,k)
    nNorm = norm(n);
    if ~isfinite(nNorm) || nNorm < 1e-12
        warning('PressureNormal:DegenerateNormal', ...
            'Degenerate pressure normal %s at load step %d substep %d.', label, k);
        n = [0 0 0];
    else
        n = n/nNorm;
    end
end







