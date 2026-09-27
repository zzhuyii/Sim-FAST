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
pressure=(32)*1000;

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
for i=1:N
% Skeleton
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
end

% Reorient and translate structural units
nodes_per_unit = 92;
move_distance  = 2*(t + t + w*cos(alpha));

% Units 1 and 2:
% Rotate -90 degrees about the global y-axis
rotation_angle = pi/2;

R_y_neg90 = [ cos(rotation_angle), 0, sin(rotation_angle);
              0,                   1, 0;
             -sin(rotation_angle), 0, cos(rotation_angle)];

for unit_id = 1:2
    node_range = (unit_id-1)*nodes_per_unit + ...
                 (1:nodes_per_unit);

    node.coordinates_mat(node_range,:) = ...
        (R_y_neg90 * node.coordinates_mat(node_range,:)')';
end


% Units 3 and 4:
% Translate in the positive x- and z-directions
for unit_id = 3:4
    node_range = (unit_id-1)*nodes_per_unit + ...
                 (1:nodes_per_unit);

    node.coordinates_mat(node_range,1) = ...
        node.coordinates_mat(node_range,1) + move_distance;

    node.coordinates_mat(node_range,3) = ...
        node.coordinates_mat(node_range,3) + move_distance;
end


% Units 5 and 6:
% Rotate +90 degrees about a line parallel to the y-axis
% Axis location:
% x = -((6-1)*(w*cos(alpha)+2*t)+t+t+w*cos(alpha))
% z = 0

rotation_angle = -pi/2;

R_y_pos90 = [ cos(rotation_angle), 0, sin(rotation_angle);
              0,                   1, 0;
             -sin(rotation_angle), 0, cos(rotation_angle)];

rotation_axis_x = ...
    -((6-1)*(w*cos(alpha)+2*t) + t + t + w*cos(alpha));

rotation_axis_point = [rotation_axis_x, 0, 0];

x_move_distance = 4 * (t + t + w*cos(alpha));

for unit_id = 5:6
    node_range = (unit_id-1)*nodes_per_unit + ...
                 (1:nodes_per_unit);

    coordinates = node.coordinates_mat(node_range,:);

    % Move the rotation axis to the origin
    coordinates = coordinates - rotation_axis_point;

    % Rotate about the y-axis
    coordinates = (R_y_pos90 * coordinates')';

    % Move the rotation axis back
    coordinates = coordinates + rotation_axis_point;

    % Translate in the positive x-direction
    coordinates(:,1) = coordinates(:,1) + x_move_distance;

    node.coordinates_mat(node_range,:) = coordinates;
end

% Frame corner cube nodes
% Frame corner cube 1
corner_nodes_1 = [
    node.coordinates_mat(181,:);        % 553
    node.coordinates_mat(182,:);        % 554
    node.coordinates_mat(183,:);        % 555
    node.coordinates_mat(184,:);        % 556
    node.coordinates_mat(187,:);        % 557
    node.coordinates_mat(188,:);        % 558
    [node.coordinates_mat(183,1), ...   
     node.coordinates_mat(183,2), ...   
     node.coordinates_mat(187,3)];      % 559
    [node.coordinates_mat(184,1), ...   
     node.coordinates_mat(184,2), ...   
     node.coordinates_mat(188,3)]       % 560
];

node.coordinates_mat = [
    node.coordinates_mat;
    corner_nodes_1
];


% Frame corner cube 2
corner_nodes_2 = [
    node.coordinates_mat(369,:);        % 561
    node.coordinates_mat(370,:);        % 562
    node.coordinates_mat(371,:);        % 563
    node.coordinates_mat(372,:);        % 564
    node.coordinates_mat(367,:);        % 565
    node.coordinates_mat(368,:);        % 566
    [node.coordinates_mat(371,1), ...  
     node.coordinates_mat(371,2), ...
     node.coordinates_mat(367,3)];      % 567
    [node.coordinates_mat(372,1), ...
     node.coordinates_mat(372,2), ...
     node.coordinates_mat(368,3)]       % 568
];

node.coordinates_mat = [
    node.coordinates_mat;
    corner_nodes_2
];

%% Copy the entire structure in the positive y-direction

f = 3;                    % Total number of structures
N=f*N;
copy_distance = 2*0.2;   % Distance between adjacent copies: 0.04 m

assert(f >= 1 && f == floor(f), ...
    'f must be a positive integer.');

% Save the original node coordinates before copying
original_coordinates = node.coordinates_mat;
nodes_per_structure = size(original_coordinates,1);

% Preallocate space for all structures
all_coordinates = zeros(f*nodes_per_structure, 3);

% Preserve the original structure
all_coordinates(1:nodes_per_structure,:) = original_coordinates;

% Generate f-1 copies
for copy_id = 2:f

    % Cumulative translation:
    % copy 1: y + 0.04
    % copy 2: y + 0.08
    % copy 3: y + 0.12
    y_translation = (copy_id-1)*copy_distance;

    copied_coordinates = original_coordinates;
    copied_coordinates(:,2) = ...
        copied_coordinates(:,2) + y_translation;

    first_node = (copy_id-1)*nodes_per_structure + 1;
    last_node  = copy_id*nodes_per_structure;

    all_coordinates(first_node:last_node,:) = copied_coordinates;
end

% Replace the node matrix with the complete copied structure
node.coordinates_mat = all_coordinates;

copy_distance_2=0.2;
% Generate Frame connections
for copy_id = 1:f-1
    y_translation = (1+(copy_id-1)*2)*copy_distance_2;
    original_cube_coordinates = node.coordinates_mat(553:568,:);
    copied_coordinates = original_cube_coordinates;
    % Move in the positive y-direction
    copied_coordinates(:,2) = copied_coordinates(:,2) + ...
                              y_translation;

    % Add the copied nodes
    node.coordinates_mat = [
        node.coordinates_mat;
        copied_coordinates
    ];
end

fprintf('Number of structures: %d\n', f);
fprintf('Nodes per structure: %d\n', nodes_per_structure);
fprintf('Total nodes: %d\n', size(node.coordinates_mat,1));

%% Define Plotting Functions
plots=Plot_Foldable_Unit;
plots.assembly=assembly;
plots.displayRange=[-0.7; 0.4; -0.5; 0.5*f; -0.1; 1];
plots.viewAngle1=20;
plots.viewAngle2=20;
plots.holdTime=0.04;

plots.Plot_Shape_Node_Number;


%% CST Define
tri_ijk=[];
tri_direction=[];
% n_e_u = size(node.coordinates_mat, 1)/N; % each unit nodes amount
n_e_u = 92;
n_e_s = 568;
n=6;
for k=1:f
for i=1:n
    % Skeleton
     tri_ijk=[tri_ijk;
        n_e_s*(k-1)+n_e_u*(i-1)+1    n_e_s*(k-1)+n_e_u*(i-1)+5    n_e_s*(k-1)+n_e_u*(i-1)+6;
        n_e_s*(k-1)+n_e_u*(i-1)+1    n_e_s*(k-1)+n_e_u*(i-1)+2    n_e_s*(k-1)+n_e_u*(i-1)+6;
        n_e_s*(k-1)+n_e_u*(i-1)+1    n_e_s*(k-1)+n_e_u*(i-1)+5    n_e_s*(k-1)+n_e_u*(i-1)+7;
        n_e_s*(k-1)+n_e_u*(i-1)+1    n_e_s*(k-1)+n_e_u*(i-1)+3    n_e_s*(k-1)+n_e_u*(i-1)+7;
        n_e_s*(k-1)+n_e_u*(i-1)+3    n_e_s*(k-1)+n_e_u*(i-1)+7    n_e_s*(k-1)+n_e_u*(i-1)+8;
        n_e_s*(k-1)+n_e_u*(i-1)+3    n_e_s*(k-1)+n_e_u*(i-1)+4    n_e_s*(k-1)+n_e_u*(i-1)+8;
        n_e_s*(k-1)+n_e_u*(i-1)+4    n_e_s*(k-1)+n_e_u*(i-1)+8    n_e_s*(k-1)+n_e_u*(i-1)+6;
        n_e_s*(k-1)+n_e_u*(i-1)+4    n_e_s*(k-1)+n_e_u*(i-1)+2    n_e_s*(k-1)+n_e_u*(i-1)+6;];  % 8

    for j=1:4
         tri_ijk=[tri_ijk;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+9    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+10    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+13;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+10    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+11    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+14;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+9    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+12    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+13;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+10    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+13    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+14;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+11    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+14    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+15;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+12    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+13    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+16;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+13    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+16    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+17;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+13    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+17    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+18;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+13    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+14    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+18;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+14    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+18    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+19;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+14    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+19    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+20;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+14    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+15    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+20;

             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+16    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+21    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+22;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+16    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+17    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+22;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+17    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+18    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+22;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+18    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+22    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+23;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+18    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+19    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+23;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+19    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+20    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+23;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+20    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+23    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+24;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+21    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+22    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+25;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+22    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+25    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+26;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+22    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+23    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+26;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+23    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+26    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+27;
             n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+23    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+24    n_e_s*(k-1)+n_e_u*(i-1)+(j-1)*19+27;];    % 104
    end
    tri_ijk=[tri_ijk;
        n_e_s*(k-1)+n_e_u*(i-1)+85    n_e_s*(k-1)+n_e_u*(i-1)+89    n_e_s*(k-1)+n_e_u*(i-1)+90;
        n_e_s*(k-1)+n_e_u*(i-1)+85    n_e_s*(k-1)+n_e_u*(i-1)+86    n_e_s*(k-1)+n_e_u*(i-1)+90;
        n_e_s*(k-1)+n_e_u*(i-1)+85    n_e_s*(k-1)+n_e_u*(i-1)+89    n_e_s*(k-1)+n_e_u*(i-1)+91;
        n_e_s*(k-1)+n_e_u*(i-1)+85    n_e_s*(k-1)+n_e_u*(i-1)+87    n_e_s*(k-1)+n_e_u*(i-1)+91;
        n_e_s*(k-1)+n_e_u*(i-1)+87    n_e_s*(k-1)+n_e_u*(i-1)+91    n_e_s*(k-1)+n_e_u*(i-1)+92;
        n_e_s*(k-1)+n_e_u*(i-1)+87    n_e_s*(k-1)+n_e_u*(i-1)+88    n_e_s*(k-1)+n_e_u*(i-1)+92;
        n_e_s*(k-1)+n_e_u*(i-1)+88    n_e_s*(k-1)+n_e_u*(i-1)+92    n_e_s*(k-1)+n_e_u*(i-1)+90;
        n_e_s*(k-1)+n_e_u*(i-1)+88    n_e_s*(k-1)+n_e_u*(i-1)+86    n_e_s*(k-1)+n_e_u*(i-1)+90;];    % 112
end

% Frame corner cube cst
tri_ijk=[tri_ijk;
    n_e_s*(k-1)+553 n_e_s*(k-1)+554 n_e_s*(k-1)+556;
    n_e_s*(k-1)+553 n_e_s*(k-1)+555 n_e_s*(k-1)+556;
    n_e_s*(k-1)+560 n_e_s*(k-1)+554 n_e_s*(k-1)+556;
    n_e_s*(k-1)+560 n_e_s*(k-1)+558 n_e_s*(k-1)+554;
    n_e_s*(k-1)+558 n_e_s*(k-1)+559 n_e_s*(k-1)+560;
    n_e_s*(k-1)+558 n_e_s*(k-1)+557 n_e_s*(k-1)+559;
    n_e_s*(k-1)+553 n_e_s*(k-1)+555 n_e_s*(k-1)+559;
    n_e_s*(k-1)+553 n_e_s*(k-1)+557 n_e_s*(k-1)+559;
    n_e_s*(k-1)+556 n_e_s*(k-1)+559 n_e_s*(k-1)+560;
    n_e_s*(k-1)+556 n_e_s*(k-1)+559 n_e_s*(k-1)+555;
    n_e_s*(k-1)+557 n_e_s*(k-1)+558 n_e_s*(k-1)+554;
    n_e_s*(k-1)+557 n_e_s*(k-1)+553 n_e_s*(k-1)+554;];

tri_ijk=[tri_ijk;
    n_e_s*(k-1)+564 n_e_s*(k-1)+562 n_e_s*(k-1)+566;
    n_e_s*(k-1)+564 n_e_s*(k-1)+566 n_e_s*(k-1)+568;
    n_e_s*(k-1)+568 n_e_s*(k-1)+566 n_e_s*(k-1)+565;
    n_e_s*(k-1)+568 n_e_s*(k-1)+567 n_e_s*(k-1)+565;
    n_e_s*(k-1)+563 n_e_s*(k-1)+561 n_e_s*(k-1)+565;
    n_e_s*(k-1)+563 n_e_s*(k-1)+567 n_e_s*(k-1)+565;
    n_e_s*(k-1)+561 n_e_s*(k-1)+562 n_e_s*(k-1)+564;
    n_e_s*(k-1)+561 n_e_s*(k-1)+563 n_e_s*(k-1)+564;
    n_e_s*(k-1)+561 n_e_s*(k-1)+562 n_e_s*(k-1)+565;
    n_e_s*(k-1)+566 n_e_s*(k-1)+562 n_e_s*(k-1)+565;
    n_e_s*(k-1)+564 n_e_s*(k-1)+567 n_e_s*(k-1)+568;
    n_e_s*(k-1)+563 n_e_s*(k-1)+564 n_e_s*(k-1)+567;];
end

% Frame connection cube cst
for l=1:f-1
    tri_ijk=[tri_ijk;
    n_e_s*f+(l-1)*16+553-552 n_e_s*f+(l-1)*16+554-552 n_e_s*f+(l-1)*16+556-552;
    n_e_s*f+(l-1)*16+553-552 n_e_s*f+(l-1)*16+555-552 n_e_s*f+(l-1)*16+556-552;
    n_e_s*f+(l-1)*16+560-552 n_e_s*f+(l-1)*16+554-552 n_e_s*f+(l-1)*16+556-552;
    n_e_s*f+(l-1)*16+560-552 n_e_s*f+(l-1)*16+558-552 n_e_s*f+(l-1)*16+554-552;
    n_e_s*f+(l-1)*16+558-552 n_e_s*f+(l-1)*16+559-552 n_e_s*f+(l-1)*16+560-552;
    n_e_s*f+(l-1)*16+558-552 n_e_s*f+(l-1)*16+557-552 n_e_s*f+(l-1)*16+559-552;
    n_e_s*f+(l-1)*16+553-552 n_e_s*f+(l-1)*16+555-552 n_e_s*f+(l-1)*16+559-552;
    n_e_s*f+(l-1)*16+553-552 n_e_s*f+(l-1)*16+557-552 n_e_s*f+(l-1)*16+559-552;
    n_e_s*f+(l-1)*16+556-552 n_e_s*f+(l-1)*16+559-552 n_e_s*f+(l-1)*16+560-552;
    n_e_s*f+(l-1)*16+556-552 n_e_s*f+(l-1)*16+559-552 n_e_s*f+(l-1)*16+555-552;
    n_e_s*f+(l-1)*16+557-552 n_e_s*f+(l-1)*16+558-552 n_e_s*f+(l-1)*16+554-552;
    n_e_s*f+(l-1)*16+557-552 n_e_s*f+(l-1)*16+553-552 n_e_s*f+(l-1)*16+554-552;];

    tri_ijk=[tri_ijk;
    n_e_s*f+(l-1)*16+564-552 n_e_s*f+(l-1)*16+562-552 n_e_s*f+(l-1)*16+566-552;
    n_e_s*f+(l-1)*16+564-552 n_e_s*f+(l-1)*16+566-552 n_e_s*f+(l-1)*16+568-552;
    n_e_s*f+(l-1)*16+568-552 n_e_s*f+(l-1)*16+566-552 n_e_s*f+(l-1)*16+565-552;
    n_e_s*f+(l-1)*16+568-552 n_e_s*f+(l-1)*16+567-552 n_e_s*f+(l-1)*16+565-552;
    n_e_s*f+(l-1)*16+563-552 n_e_s*f+(l-1)*16+561-552 n_e_s*f+(l-1)*16+565-552;
    n_e_s*f+(l-1)*16+563-552 n_e_s*f+(l-1)*16+567-552 n_e_s*f+(l-1)*16+565-552;
    n_e_s*f+(l-1)*16+561-552 n_e_s*f+(l-1)*16+562-552 n_e_s*f+(l-1)*16+564-552;
    n_e_s*f+(l-1)*16+561-552 n_e_s*f+(l-1)*16+563-552 n_e_s*f+(l-1)*16+564-552;
    n_e_s*f+(l-1)*16+561-552 n_e_s*f+(l-1)*16+562-552 n_e_s*f+(l-1)*16+565-552;
    n_e_s*f+(l-1)*16+566-552 n_e_s*f+(l-1)*16+562-552 n_e_s*f+(l-1)*16+565-552;
    n_e_s*f+(l-1)*16+564-552 n_e_s*f+(l-1)*16+567-552 n_e_s*f+(l-1)*16+568-552;
    n_e_s*f+(l-1)*16+563-552 n_e_s*f+(l-1)*16+564-552 n_e_s*f+(l-1)*16+567-552;];
end

cst.node_ijk_mat=tri_ijk;
cstNum=size(cst.node_ijk_mat,1);

cst.t_vec=0.006*ones(cstNum,1);
cst.E_vec=4*10^9*ones(cstNum,1);
cst.v_vec=0.2*ones(cstNum,1);

plots.Plot_Shape_CST_Number;    



%% Define Rotational Spring
for k=1:f
for i=1:n
    % Skeleton
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        n_e_s*(k-1)+n_e_u*(i-1)+6    n_e_s*(k-1)+n_e_u*(i-1)+1    n_e_s*(k-1)+n_e_u*(i-1)+5    n_e_s*(k-1)+n_e_u*(i-1)+7;
        n_e_s*(k-1)+n_e_u*(i-1)+1    n_e_s*(k-1)+n_e_u*(i-1)+3    n_e_s*(k-1)+n_e_u*(i-1)+7    n_e_s*(k-1)+n_e_u*(i-1)+8;
        n_e_s*(k-1)+n_e_u*(i-1)+3    n_e_s*(k-1)+n_e_u*(i-1)+4    n_e_s*(k-1)+n_e_u*(i-1)+8    n_e_s*(k-1)+n_e_u*(i-1)+6;
        n_e_s*(k-1)+n_e_u*(i-1)+4    n_e_s*(k-1)+n_e_u*(i-1)+2    n_e_s*(k-1)+n_e_u*(i-1)+6    n_e_s*(k-1)+n_e_u*(i-1)+1;
        n_e_s*(k-1)+n_e_u*(i-1)+5    n_e_s*(k-1)+n_e_u*(i-1)+1    n_e_s*(k-1)+n_e_u*(i-1)+7    n_e_s*(k-1)+n_e_u*(i-1)+3;
        n_e_s*(k-1)+n_e_u*(i-1)+7    n_e_s*(k-1)+n_e_u*(i-1)+3    n_e_s*(k-1)+n_e_u*(i-1)+8    n_e_s*(k-1)+n_e_u*(i-1)+4;
        n_e_s*(k-1)+n_e_u*(i-1)+8    n_e_s*(k-1)+n_e_u*(i-1)+4    n_e_s*(k-1)+n_e_u*(i-1)+6    n_e_s*(k-1)+n_e_u*(i-1)+2;
        n_e_s*(k-1)+n_e_u*(i-1)+5    n_e_s*(k-1)+n_e_u*(i-1)+1    n_e_s*(k-1)+n_e_u*(i-1)+6    n_e_s*(k-1)+n_e_u*(i-1)+2;];

        for j=1:4
        rot_spr_4N.node_ijkl_mat=[
            rot_spr_4N.node_ijkl_mat;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+9    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+10    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+10    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+11;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+10    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+11    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+15;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+11    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+15    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+20;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+15    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+20    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+19;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+20    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+19    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+19    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+10    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+17;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+17    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+16;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+17    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+16    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+12;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+16    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+12    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+9;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+10    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+9    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+12;

            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+16    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+17    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+17    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+19;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+19    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+20;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+19    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+20    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+24;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+20    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+24    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+27;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+24    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+27    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+26;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+27    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+26    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+26;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+26    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+25;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+26    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+25    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+21;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+25    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+21    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+16;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+21    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+16    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+17;];
        end
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        n_e_s*(k-1)+n_e_u*(i-1)+90    n_e_s*(k-1)+n_e_u*(i-1)+85    n_e_s*(k-1)+n_e_u*(i-1)+89    n_e_s*(k-1)+n_e_u*(i-1)+91;
        n_e_s*(k-1)+n_e_u*(i-1)+85    n_e_s*(k-1)+n_e_u*(i-1)+87    n_e_s*(k-1)+n_e_u*(i-1)+91    n_e_s*(k-1)+n_e_u*(i-1)+92;
        n_e_s*(k-1)+n_e_u*(i-1)+87    n_e_s*(k-1)+n_e_u*(i-1)+88    n_e_s*(k-1)+n_e_u*(i-1)+92    n_e_s*(k-1)+n_e_u*(i-1)+90;
        n_e_s*(k-1)+n_e_u*(i-1)+88    n_e_s*(k-1)+n_e_u*(i-1)+86    n_e_s*(k-1)+n_e_u*(i-1)+90    n_e_s*(k-1)+n_e_u*(i-1)+85;
        n_e_s*(k-1)+n_e_u*(i-1)+89    n_e_s*(k-1)+n_e_u*(i-1)+85    n_e_s*(k-1)+n_e_u*(i-1)+91    n_e_s*(k-1)+n_e_u*(i-1)+87;
        n_e_s*(k-1)+n_e_u*(i-1)+91    n_e_s*(k-1)+n_e_u*(i-1)+87    n_e_s*(k-1)+n_e_u*(i-1)+92    n_e_s*(k-1)+n_e_u*(i-1)+88;
        n_e_s*(k-1)+n_e_u*(i-1)+92    n_e_s*(k-1)+n_e_u*(i-1)+88    n_e_s*(k-1)+n_e_u*(i-1)+90    n_e_s*(k-1)+n_e_u*(i-1)+86;
        n_e_s*(k-1)+n_e_u*(i-1)+86    n_e_s*(k-1)+n_e_u*(i-1)+85    n_e_s*(k-1)+n_e_u*(i-1)+90    n_e_s*(k-1)+n_e_u*(i-1)+89;]; % 120

        for j=1:4
        rot_spr_4N_D.node_ijkl_mat=[
            rot_spr_4N_D.node_ijkl_mat;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+16    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+17    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+13    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+17    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+22;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+18    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+19    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23;
            n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+14    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+19    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+20    n_e_s*(k-1)+n_e_u*(i-1)+19*(j-1)+23;]; % 136     Folding line 
        end
end

rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        n_e_s*(k-1)+555  n_e_s*(k-1)+553  n_e_s*(k-1)+556  n_e_s*(k-1)+554;
        n_e_s*(k-1)+556  n_e_s*(k-1)+560  n_e_s*(k-1)+554  n_e_s*(k-1)+558;
        n_e_s*(k-1)+560  n_e_s*(k-1)+558  n_e_s*(k-1)+559  n_e_s*(k-1)+557;
        n_e_s*(k-1)+557  n_e_s*(k-1)+559  n_e_s*(k-1)+553  n_e_s*(k-1)+555;
        n_e_s*(k-1)+555  n_e_s*(k-1)+556  n_e_s*(k-1)+559  n_e_s*(k-1)+560;
        n_e_s*(k-1)+553  n_e_s*(k-1)+554  n_e_s*(k-1)+557  n_e_s*(k-1)+558;
        n_e_s*(k-1)+553  n_e_s*(k-1)+554  n_e_s*(k-1)+556  n_e_s*(k-1)+560;
        n_e_s*(k-1)+554  n_e_s*(k-1)+558  n_e_s*(k-1)+560  n_e_s*(k-1)+559;
        n_e_s*(k-1)+558  n_e_s*(k-1)+557  n_e_s*(k-1)+559  n_e_s*(k-1)++553;
        n_e_s*(k-1)+559  n_e_s*(k-1)+555  n_e_s*(k-1)+553  n_e_s*(k-1)+556;
        n_e_s*(k-1)+554  n_e_s*(k-1)+556  n_e_s*(k-1)+560  n_e_s*(k-1)+559;
        n_e_s*(k-1)+560  n_e_s*(k-1)+554  n_e_s*(k-1)+558  n_e_s*(k-1)+557;
        n_e_s*(k-1)+554  n_e_s*(k-1)+557  n_e_s*(k-1)+553  n_e_s*(k-1)+559;
        n_e_s*(k-1)+553  n_e_s*(k-1)+559  n_e_s*(k-1)+555  n_e_s*(k-1)+556;
        n_e_s*(k-1)+553  n_e_s*(k-1)+555  n_e_s*(k-1)+556  n_e_s*(k-1)+559;
        n_e_s*(k-1)+556  n_e_s*(k-1)+559  n_e_s*(k-1)+560  n_e_s*(k-1)+558;
        n_e_s*(k-1)+559  n_e_s*(k-1)+557  n_e_s*(k-1)+558  n_e_s*(k-1)+554;
        n_e_s*(k-1)+557  n_e_s*(k-1)+553  n_e_s*(k-1)+554  n_e_s*(k-1)+556;

        n_e_s*(k-1)+563  n_e_s*(k-1)+561  n_e_s*(k-1)+564  n_e_s*(k-1)+562;
        n_e_s*(k-1)+568  n_e_s*(k-1)+566  n_e_s*(k-1)+564  n_e_s*(k-1)+562;
        n_e_s*(k-1)+567  n_e_s*(k-1)+565  n_e_s*(k-1)+568  n_e_s*(k-1)+566;
        n_e_s*(k-1)+567  n_e_s*(k-1)+565  n_e_s*(k-1)+563  n_e_s*(k-1)+561;
        n_e_s*(k-1)+561  n_e_s*(k-1)+562  n_e_s*(k-1)+565  n_e_s*(k-1)+566;
        n_e_s*(k-1)+563  n_e_s*(k-1)+564  n_e_s*(k-1)+567  n_e_s*(k-1)+568;
        n_e_s*(k-1)+561  n_e_s*(k-1)+562  n_e_s*(k-1)+564  n_e_s*(k-1)+566;
        n_e_s*(k-1)+564  n_e_s*(k-1)+566  n_e_s*(k-1)+568  n_e_s*(k-1)+565;
        n_e_s*(k-1)+568  n_e_s*(k-1)+565  n_e_s*(k-1)+567  n_e_s*(k-1)+563;
        n_e_s*(k-1)+565  n_e_s*(k-1)+561  n_e_s*(k-1)+563  n_e_s*(k-1)+564;
        n_e_s*(k-1)+565  n_e_s*(k-1)+566  n_e_s*(k-1)+562  n_e_s*(k-1)+564;
        n_e_s*(k-1)+566  n_e_s*(k-1)+564  n_e_s*(k-1)+568  n_e_s*(k-1)+567;
        n_e_s*(k-1)+564  n_e_s*(k-1)+567  n_e_s*(k-1)+563  n_e_s*(k-1)+565;
        n_e_s*(k-1)+563  n_e_s*(k-1)+565  n_e_s*(k-1)+561  n_e_s*(k-1)+562;
        n_e_s*(k-1)+564  n_e_s*(k-1)+561  n_e_s*(k-1)+562  n_e_s*(k-1)+565;
        n_e_s*(k-1)+562  n_e_s*(k-1)+565  n_e_s*(k-1)+566  n_e_s*(k-1)+568;
        n_e_s*(k-1)+565  n_e_s*(k-1)+561  n_e_s*(k-1)+563  n_e_s*(k-1)+564;
        n_e_s*(k-1)+565  n_e_s*(k-1)+567  n_e_s*(k-1)+568  n_e_s*(k-1)+564;
        n_e_s*(k-1)+567  n_e_s*(k-1)+563  n_e_s*(k-1)+564  n_e_s*(k-1)+561;];

end

% Frame connection cube 4-node-spring
for l=1:f-1
rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        n_e_s*f+(l-1)*16+555-552  n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+556-552  n_e_s*f+(l-1)*16+554-552;
        n_e_s*f+(l-1)*16+556-552  n_e_s*f+(l-1)*16+560-552  n_e_s*f+(l-1)*16+554-552  n_e_s*f+(l-1)*16+558-552;
        n_e_s*f+(l-1)*16+560-552  n_e_s*f+(l-1)*16+558-552  n_e_s*f+(l-1)*16+559-552  n_e_s*f+(l-1)*16+557-552;
        n_e_s*f+(l-1)*16+557-552  n_e_s*f+(l-1)*16+559-552  n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+555-552;
        n_e_s*f+(l-1)*16+555-552  n_e_s*f+(l-1)*16+556-552  n_e_s*f+(l-1)*16+559-552  n_e_s*f+(l-1)*16+560-552;
        n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+554-552  n_e_s*f+(l-1)*16+557-552  n_e_s*f+(l-1)*16+558-552;
        n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+554-552  n_e_s*f+(l-1)*16+556-552  n_e_s*f+(l-1)*16+560-552;
        n_e_s*f+(l-1)*16+554-552  n_e_s*f+(l-1)*16+558-552  n_e_s*f+(l-1)*16+560-552  n_e_s*f+(l-1)*16+559-552;
        n_e_s*f+(l-1)*16+558-552  n_e_s*f+(l-1)*16+557-552  n_e_s*f+(l-1)*16+559-552  n_e_s*f+(l-1)*16+553-552;
        n_e_s*f+(l-1)*16+559-552  n_e_s*f+(l-1)*16+555-552  n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+556-552;
        n_e_s*f+(l-1)*16+554-552  n_e_s*f+(l-1)*16+556-552  n_e_s*f+(l-1)*16+560-552  n_e_s*f+(l-1)*16+559-552;
        n_e_s*f+(l-1)*16+560-552  n_e_s*f+(l-1)*16+554-552  n_e_s*f+(l-1)*16+558-552  n_e_s*f+(l-1)*16+557-552;
        n_e_s*f+(l-1)*16+554-552  n_e_s*f+(l-1)*16+557-552  n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+559-552;
        n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+559-552  n_e_s*f+(l-1)*16+555-552  n_e_s*f+(l-1)*16+556-552;
        n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+555-552  n_e_s*f+(l-1)*16+556-552  n_e_s*f+(l-1)*16+559-552;
        n_e_s*f+(l-1)*16+556-552  n_e_s*f+(l-1)*16+559-552  n_e_s*f+(l-1)*16+560-552  n_e_s*f+(l-1)*16+558-552;
        n_e_s*f+(l-1)*16+559-552  n_e_s*f+(l-1)*16+557-552  n_e_s*f+(l-1)*16+558-552  n_e_s*f+(l-1)*16+554-552;
        n_e_s*f+(l-1)*16+557-552  n_e_s*f+(l-1)*16+553-552  n_e_s*f+(l-1)*16+554-552  n_e_s*f+(l-1)*16+556-552;

        n_e_s*f+(l-1)*16+563-552  n_e_s*f+(l-1)*16+561-552  n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+562-552;
        n_e_s*f+(l-1)*16+568-552  n_e_s*f+(l-1)*16+566-552  n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+562-552;
        n_e_s*f+(l-1)*16+567-552  n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+568-552  n_e_s*f+(l-1)*16+566-552;
        n_e_s*f+(l-1)*16+567-552  n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+563-552  n_e_s*f+(l-1)*16+561-552;
        n_e_s*f+(l-1)*16+561-552  n_e_s*f+(l-1)*16+562-552  n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+566-552;
        n_e_s*f+(l-1)*16+563-552  n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+567-552  n_e_s*f+(l-1)*16+568-552;
        n_e_s*f+(l-1)*16+561-552  n_e_s*f+(l-1)*16+562-552  n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+566-552;
        n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+566-552  n_e_s*f+(l-1)*16+568-552  n_e_s*f+(l-1)*16+565-552;
        n_e_s*f+(l-1)*16+568-552  n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+567-552  n_e_s*f+(l-1)*16+563-552;
        n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+561-552  n_e_s*f+(l-1)*16+563-552  n_e_s*f+(l-1)*16+564-552;
        n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+566-552  n_e_s*f+(l-1)*16+562-552  n_e_s*f+(l-1)*16+564-552;
        n_e_s*f+(l-1)*16+566-552  n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+568-552  n_e_s*f+(l-1)*16+567-552;
        n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+567-552  n_e_s*f+(l-1)*16+563-552  n_e_s*f+(l-1)*16+565-552;
        n_e_s*f+(l-1)*16+563-552  n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+561-552  n_e_s*f+(l-1)*16+562-552;
        n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+561-552  n_e_s*f+(l-1)*16+562-552  n_e_s*f+(l-1)*16+565-552;
        n_e_s*f+(l-1)*16+562-552  n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+566-552  n_e_s*f+(l-1)*16+568-552;
        n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+561-552  n_e_s*f+(l-1)*16+563-552  n_e_s*f+(l-1)*16+564-552;
        n_e_s*f+(l-1)*16+565-552  n_e_s*f+(l-1)*16+567-552  n_e_s*f+(l-1)*16+568-552  n_e_s*f+(l-1)*16+564-552;
        n_e_s*f+(l-1)*16+567-552  n_e_s*f+(l-1)*16+563-552  n_e_s*f+(l-1)*16+564-552  n_e_s*f+(l-1)*16+561-552;];
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
rot_spr_4N_D.mv_factor_vec = 100 * ones(rot_dir_num, 1);

% The Mountain Valley Assignment for the directional spring
rot_spr_4N_D.mv_vec = ones(rot_dir_num,1);

for j=1:f
for i = 1:n
    idx = (i-1)*16 + (j-1)*16*6 + [1:4, 9:12];   
    rot_spr_4N_D.mv_vec(idx) = 0;
end
end

plots.Plot_Shape_Spr_Number();
plots.Plot_Shape_DirectSpr_Number();


%% Define the connectors
zlsprStiff=10000000;
for k=1:f

for i=1:n 
    zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+(i-1)*n_e_u+5  n_e_s*(k-1)+(i-1)*n_e_u+9;  n_e_s*(k-1)+(i-1)*n_e_u+5  n_e_s*(k-1)+(i-1)*n_e_u+28];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+(i-1)*n_e_u+7  n_e_s*(k-1)+(i-1)*n_e_u+30;  n_e_s*(k-1)+(i-1)*n_e_u+7  n_e_s*(k-1)+(i-1)*n_e_u+47];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+(i-1)*n_e_u+8  n_e_s*(k-1)+(i-1)*n_e_u+49;  n_e_s*(k-1)+(i-1)*n_e_u+8  n_e_s*(k-1)+(i-1)*n_e_u+68];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+(i-1)*n_e_u+6  n_e_s*(k-1)+(i-1)*n_e_u+66;  n_e_s*(k-1)+(i-1)*n_e_u+6  n_e_s*(k-1)+(i-1)*n_e_u+11];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+(i-1)*n_e_u+85  n_e_s*(k-1)+(i-1)*n_e_u+25;  n_e_s*(k-1)+(i-1)*n_e_u+85  n_e_s*(k-1)+(i-1)*n_e_u+44];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+(i-1)*n_e_u+87  n_e_s*(k-1)+(i-1)*n_e_u+46;  n_e_s*(k-1)+(i-1)*n_e_u+87  n_e_s*(k-1)+(i-1)*n_e_u+63];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+(i-1)*n_e_u+88  n_e_s*(k-1)+(i-1)*n_e_u+65;  n_e_s*(k-1)+(i-1)*n_e_u+88  n_e_s*(k-1)+(i-1)*n_e_u+84];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+(i-1)*n_e_u+86  n_e_s*(k-1)+(i-1)*n_e_u+82;  n_e_s*(k-1)+(i-1)*n_e_u+86  n_e_s*(k-1)+(i-1)*n_e_u+27];
end

zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+89  n_e_s*(k-1)+93];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+90  n_e_s*(k-1)+94];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+92  n_e_s*(k-1)+96];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+91  n_e_s*(k-1)+95];

zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+181  n_e_s*(k-1)+553];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+182  n_e_s*(k-1)+554];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+183  n_e_s*(k-1)+555];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+184  n_e_s*(k-1)+556];

zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+553  n_e_s*(k-1)+185];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+554  n_e_s*(k-1)+186];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+557  n_e_s*(k-1)+187];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+558  n_e_s*(k-1)+188];

zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+273  n_e_s*(k-1)+277];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+274  n_e_s*(k-1)+278];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+275  n_e_s*(k-1)+279];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+276  n_e_s*(k-1)+280];

zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+365  n_e_s*(k-1)+561];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+366  n_e_s*(k-1)+562];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+367  n_e_s*(k-1)+565];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+368  n_e_s*(k-1)+566];

zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+369  n_e_s*(k-1)+561];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+370  n_e_s*(k-1)+562];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+371  n_e_s*(k-1)+563];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+372  n_e_s*(k-1)+564];

zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+461  n_e_s*(k-1)+457];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+462  n_e_s*(k-1)+458];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+463  n_e_s*(k-1)+459];
zlspr.node_ij_mat=[zlspr.node_ij_mat;n_e_s*(k-1)+464  n_e_s*(k-1)+460];

end

for j=1:f-1
    zlspr.node_ij_mat=[zlspr.node_ij_mat;553+(j-1)*568  n_e_s*f+554-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;555+(j-1)*568  n_e_s*f+556-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;557+(j-1)*568  n_e_s*f+558-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;559+(j-1)*568  n_e_s*f+560-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;561+(j-1)*568  n_e_s*f+562-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;563+(j-1)*568  n_e_s*f+564-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;565+(j-1)*568  n_e_s*f+566-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;567+(j-1)*568  n_e_s*f+568-552+16*(j-1)];

    zlspr.node_ij_mat=[zlspr.node_ij_mat;554+(j)*568  n_e_s*f+553-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;556+(j)*568  n_e_s*f+555-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;558+(j)*568  n_e_s*f+557-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;560+(j)*568  n_e_s*f+559-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;562+(j)*568  n_e_s*f+561-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;564+(j)*568  n_e_s*f+563-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;566+(j)*568  n_e_s*f+565-552+16*(j-1)];
    zlspr.node_ij_mat=[zlspr.node_ij_mat;568+(j)*568  n_e_s*f+567-552+16*(j-1)];
end
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
    1 1 1 1;    
    2 1 1 1;
    3 1 1 1;
    4 1 1 1;
    549 0 1 1;
    550 0 1 1;
    551 0 1 1;
    552 0 1 1;];

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
    equal_nodal_force_direct = zeros(8*n*f,4);
    equal_nodal_force_direct(:,1) = (1:8*n*f)';
    for l=1:f
    for i=1:n
        
        n1=cross((xcurrent(19+n_e_u*(i-1)+568*(l-1),:)-xcurrent(17+n_e_u*(i-1)+568*(l-1),:)),(xcurrent(13+n_e_u*(i-1)+568*(l-1),:)-xcurrent(17+n_e_u*(i-1)+568*(l-1),:)));
        n2=cross((xcurrent(22+n_e_u*(i-1)+568*(l-1),:)-xcurrent(17+n_e_u*(i-1)+568*(l-1),:)),(xcurrent(19+n_e_u*(i-1)+568*(l-1),:)-xcurrent(17+n_e_u*(i-1)+568*(l-1),:)));

        n3=cross((xcurrent(38+n_e_u*(i-1)+568*(l-1),:)-xcurrent(36+n_e_u*(i-1)+568*(l-1),:)),(xcurrent(32+n_e_u*(i-1)+568*(l-1),:)-xcurrent(36+n_e_u*(i-1)+568*(l-1),:)));
        n4=cross((xcurrent(41+n_e_u*(i-1)+568*(l-1),:)-xcurrent(36+n_e_u*(i-1)+568*(l-1),:)),(xcurrent(38+n_e_u*(i-1)+568*(l-1),:)-xcurrent(36+n_e_u*(i-1)+568*(l-1),:)));

        n5=cross((xcurrent(51+n_e_u*(i-1)+568*(l-1),:)-xcurrent(55+n_e_u*(i-1)+568*(l-1),:)),(xcurrent(57+n_e_u*(i-1)+568*(l-1),:)-xcurrent(55+n_e_u*(i-1)+568*(l-1),:)));
        n6=cross((xcurrent(57+n_e_u*(i-1)+568*(l-1),:)-xcurrent(55+n_e_u*(i-1)+568*(l-1),:)),(xcurrent(60+n_e_u*(i-1)+568*(l-1),:)-xcurrent(55+n_e_u*(i-1)+568*(l-1),:)));

        n7=cross((xcurrent(70+n_e_u*(i-1)+568*(l-1),:)-xcurrent(74+n_e_u*(i-1)+568*(l-1),:)),(xcurrent(76+n_e_u*(i-1)+568*(l-1),:)-xcurrent(74+n_e_u*(i-1)+568*(l-1),:)));
        n8=cross((xcurrent(76+n_e_u*(i-1)+568*(l-1),:)-xcurrent(74+n_e_u*(i-1)+568*(l-1),:)),(xcurrent(79+n_e_u*(i-1)+568*(l-1),:)-xcurrent(74+n_e_u*(i-1)+568*(l-1),:)));

        n1=safeUnitNormal(n1,'n1',k);
        n2=safeUnitNormal(n2,'n2',k);
        n3=safeUnitNormal(n3,'n3',k);
        n4=safeUnitNormal(n4,'n4',k);
        n5=safeUnitNormal(n5,'n5',k);
        n6=safeUnitNormal(n6,'n6',k);
        n7=safeUnitNormal(n7,'n7',k);
        n8=safeUnitNormal(n8,'n8',k);

        equal_nodal_force_direct(1+8*(i-1)+8*n*(l-1),2:4)=n1;
        equal_nodal_force_direct(2+8*(i-1)+8*n*(l-1),2:4)=n2;
        equal_nodal_force_direct(3+8*(i-1)+8*n*(l-1),2:4)=n3;
        equal_nodal_force_direct(4+8*(i-1)+8*n*(l-1),2:4)=n4;
        equal_nodal_force_direct(5+8*(i-1)+8*n*(l-1),2:4)=n5;
        equal_nodal_force_direct(6+8*(i-1)+8*n*(l-1),2:4)=n6;
        equal_nodal_force_direct(7+8*(i-1)+8*n*(l-1),2:4)=n7;
        equal_nodal_force_direct(8+8*(i-1)+8*n*(l-1),2:4)=n8;
    end

    % Nodal force apply
    nodeNum=size(node.coordinates_mat,1);
    nr.load = zeros(nodeNum, 4);
    nr.load(:,1) = (1:nodeNum)';
    for l=1:f
    for i=1:n
        for j=1:4
        nr.load(13+19*(j-1)+n_e_u*(i-1)+568*(l-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*equal_nodal_force_direct(1+2*(j-1)+8*(i-1)+8*n*(l-1),2:4);
        nr.load(14+19*(j-1)+n_e_u*(i-1)+568*(l-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*equal_nodal_force_direct(1+2*(j-1)+8*(i-1)+8*n*(l-1),2:4);
        nr.load(17+19*(j-1)+n_e_u*(i-1)+568*(l-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*-equal_nodal_force_direct(1+2*(j-1)+8*(i-1)+8*n*(l-1),2:4)+((k)*(equal_nodal_force_magnitude/step))*-equal_nodal_force_direct(2+2*(j-1)+8*(i-1)+8*n*(l-1),2:4);
        nr.load(19+19*(j-1)+n_e_u*(i-1)+568*(l-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*-equal_nodal_force_direct(1+2*(j-1)+8*(i-1)+8*n*(l-1),2:4)+((k)*(equal_nodal_force_magnitude/step))*-equal_nodal_force_direct(2+2*(j-1)+8*(i-1)+8*n*(l-1),2:4);
        nr.load(22+19*(j-1)+n_e_u*(i-1)+568*(l-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*equal_nodal_force_direct(2+2*(j-1)+8*(i-1)+8*n*(l-1),2:4);
        nr.load(23+19*(j-1)+n_e_u*(i-1)+568*(l-1), 2:4) = ((k)*(equal_nodal_force_magnitude/step))*equal_nodal_force_direct(2+2*(j-1)+8*(i-1)+8*n*(l-1),2:4);
        end
    end
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

plots.fileName='Big_Frame1.gif';
plots.Plot_Deformed_His(Uhis(1:2:end,:,:))


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







