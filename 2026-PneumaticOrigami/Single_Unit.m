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
l_x=0.04; 
l_y=0.08;

% Layers amount
layer_num = 6; 

% each layer's angle
square_angle = (1/layer_num) * (pi-2*alpha); 

% Unit Number
N=1; 

% Target pressure to be applied
pressure=(21.1)*1000;

% Target load on the structure
External_load_mass=0; % kg
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


%% Define Plotting Functions
plots=Plot_Foldable_Unit;
plots.assembly=assembly;
plots.displayRange=[-1*w*(N+1); 0.1; -0.1; 0.3; -0.1; 0.3];
plots.viewAngle1=20;
plots.viewAngle2=20;
plots.holdTime=0.04;

plots.Plot_Shape_Node_Number;


%% CST Define
tri_ijk=[];
tri_direction=[];
n_e_u = size(node.coordinates_mat, 1)/N; % each unit nodes amount

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
rot_spr_4N_D.mv_factor_vec=100*ones(rot_dir_num,1);

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
if N<=1
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
else
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

    for i=1:N-1 
        zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+89 (i)*n_e_u+1];
        zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+90 (i)*n_e_u+2];
        zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+91 (i)*n_e_u+3];
        zlspr.node_ij_mat=[zlspr.node_ij_mat;(i-1)*n_e_u+92 (i)*n_e_u+4];
    end
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
    89 0 0 1;
    90 0 0 1;
    91 0 0 1;
    92 0 0 1;];

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
    
    % load put on the structure
    nr.load(89, 2:4) = nr.load(89, 2:4) + [(k)*(load_on_structure/step) 0 0];
    nr.load(90, 2:4) = nr.load(90, 2:4) + [(k)*(load_on_structure/step) 0 0];
    nr.load(91, 2:4) = nr.load(91, 2:4) + [(k)*(load_on_structure/step) 0 0];
    nr.load(92, 2:4) = nr.load(92, 2:4) + [(k)*(load_on_structure/step) 0 0];    

    Uhis(k,:,:)=squeeze(nr.Solve());

end

plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)))

plots.fileName='Single_Unit.gif';
plots.Plot_Deformed_His(Uhis(1:2:end,:,:))


%% Angle measurement
p31 = xcurrent(31,:);
p35 = xcurrent(35,:);
p36 = xcurrent(36,:);
p40 = xcurrent(40,:);
a_m_n_1 = cross(p35 - p31, p36 - p31);
a_m_n_2 = cross(p36 - p35, p40 - p35);
cos_theta = dot(a_m_n_1, a_m_n_2) / (norm(a_m_n_1) * norm(a_m_n_2));
cos_theta = max(-1,min(1,cos_theta));
theta_rad = acos(cos_theta);
theta_deg = rad2deg(theta_rad); 


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







