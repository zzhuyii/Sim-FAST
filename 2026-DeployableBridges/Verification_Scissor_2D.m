clear all;
clc;
close all;

% Geometry of the two unit deployable scissor bridge in citation:
% Y. Chikahiro, I. Ario, M. Nakazawa, S. Ono, J. Holnicki-Szulc,
% P. Pawlowski, C. Graczykowski, A. Watson, Experimental and numerical
% study of full-scale scissor type bridge, 
% Automation in Construction 71 (2016) 171–180.

H=2; % Height (m) 
L=3.5; % Length (m)

% The members of the main frame
barA_A6N01=0.0028; % A: Cross section area (m^2)
barE_A6N01=6.25*10^10; % E: Youngs' modulus (Pa)
barI_A6N01=1.1463*10^-5; % I: Moment of Inertia (m^4) 
barSigmaY_A6N01=1.8*10^8; % σy: Yield Strength (Pa)

% The deck members
deckA_A6063=0.00304; 
deckE_A6063=6.8*10^10; 
deckI_A6063=1.827*10^-6; 
deckSigmaY_A6063=1.1*10^8;

% member information
barL=sqrt(H^2+L^2);
kspr_bar=3*barE_A6N01*barI_A6N01/barL;
kspr_deck=3*deckE_A6063*deckI_A6063/(2*L);

% Define all the objects for the elements we need to use
node=Elements_Nodes;
bar=CD_Elements_Bars;
rot_spr_3N=CD_Elements_RotSprings_3N;

% Create tje truss Assembly
assembly=Assembly_2D_Mechanism;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_3N=rot_spr_3N;


% Loading cases
PF=5.2e3;  % N  Front axle
PR=4.4e3;  % N  Rear axle
s=3.1; 

xR=L-s/2; % Rear axle's position
xF=L+s/2; % Front axle's position


%% Define the nodal coordinates
% The Y-direction coordinates are all zero.
node.coordinates_mat=[
    0      0   0;
    L      0   0;
    2*L    0   0;
    0      0   H;
    L      0   H;
    2*L    0   H;
    L/2    0   H/2;
    3*L/2  0   H/2;
    L/4    0   H/4;
    3*L/4  0   H/4;
    5*L/4  0   H/4;
    7*L/4  0   H/4;
    L/4    0   3*H/4;
    3*L/4  0   3*H/4;
    5*L/4  0   3*H/4;
    7*L/4  0   3*H/4;
    xR     0   0; 
    xF     0   0; 
    ];

nR=size(node.coordinates_mat,1);
nF=nR-1;

% Set up the plotting function for inspection
plots=Plot_2D_Mechanism();
plots.assembly=assembly;

% Range of the plot
plots.displayRange=[-1;2*L+1;-1;1;-1;H+0.5]; 

% Change the viewing angle
plots.viewAngle1=0;
plots.viewAngle2=0;

% Plot the nodal coordinates for inspection
plots.Plot_Shape_Node_Number()


%% Define how bars are connected
bar.node_ij_mat=[
    1    9;
    9    7;
    7    14;
    14   5;
    4    13;
    13   7;
    7    10;
    10   2;
    
    2    11;
    11   8;
    8    16;
    16   6;
    5    15;
    15   8;
    8    12;
    12   3;
    ];

% Define the bottom deck % ******
bar.node_ij_mat=[bar.node_ij_mat;
    1   17;
    17  2;
    2   18;
    18  3;
    ];

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA_A6N01*ones(barNum,1);
bar.E_vec=barE_A6N01*ones(barNum,1);

bar.A_vec=[barA_A6N01*ones(barNum-6,1);  deckA_A6063*ones(6,1)]; % ******
bar.E_vec=[barE_A6N01*ones(barNum-6,1);  deckE_A6063*ones(6,1)]; % ******
 
plots.Plot_Shape_Bar_Number();


%% Define how rotational springs are connected
rot_spr_3N.node_ijk_mat=[
    1   9   7;
    9   7   14;
    7   14  5;

    4   13  7;
    13  7   10;
    7   10  2;

    2   11  8;
    11  8   16;
    8   16  6;

    5   15  8;
    15  8   12;
    8   12  3; % 12
    ];

% Define the bottom deck 
rot_spr_3N.node_ijk_mat = [rot_spr_3N.node_ijk_mat;
    1   17  2;
    2  18  3; 
    ];

rot_spr_3N.rot_spr_K_vec = [kspr_bar*ones(12,1); kspr_deck*ones(2,1)]; 

% Initialize the entire assembly again for the new rot-spring
assembly.Initialize_Assembly();
plots.Plot_Shape_Spr_Number();


%% Set up the loading solver
nr=Solver_NR_Loading;
nr.assembly=assembly;
nr.supp=[1,1,1,1; % Supported
         2,0,1,0;
         3,0,1,1; % Supported
         4,0,1,0;
         5,0,1,0;
         6,0,1,0;
         7,0,1,0;
         8,0,1,0;
         9,0,1,0;
         10,0,1,0;
         11,0,1,0;
         12,0,1,0;
         13,0,1,0;
         14,0,1,0;
         15,0,1,0;
         16,0,1,0;
         17,0,1,0;
         18,0,1,0; 
         ];

% Set up the load 
% bridge self weight 9600N 
% vehicle load 11800N, 13800N
% load is distributed to two scissors
% self weight uniformly distributed 

nr.load=[1,0,0,-9600/2/18/5;
         2,0,0,-9600/2/18/5;
         3,0,0,-9600/2/18/5;
         4,0,0,-9600/2/18/5;
         5,0,0,-9600/2/18/5;
         6,0,0,-9600/2/18/5;
         7,0,0,-9600/2/18/5;
         8,0,0,-9600/2/18/5;
         9,0,0,-9600/2/18/5;
         10,0,0,-9600/2/18/5;
         11,0,0,-9600/2/18/5;
         12,0,0,-9600/2/18/5;
         13,0,0,-9600/2/18/5;
         14,0,0,-9600/2/18/5;
         15,0,0,-9600/2/18/5;
         16,0,0,-9600/2/18/5;
         17,0,0,-4400/2/5-9600/2/18/5; % R
         18,0,0,-5200/2/5-9600/2/18/5; % F
         ]; 

% Set up the total loading step
nr.increStep=5;
% Set up the maximum iteration
nr.iterMax=30;
% Set up the tolorence
nr.tol=10^-6;

% Solve for the deformation history
Uhis=nr.Solve();

% Plot the deformed shape
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)),zeros(size(squeeze(Uhis(end,:,:)))));

% Find the moment of the rotational springs
rot_spr_3N.Solve_Global_Theta(node,squeeze(Uhis(end,:,:)));
moment_vec = rot_spr_3N.rot_spr_K_vec.*(rot_spr_3N.theta_current_vec-rot_spr_3N.theta_stress_free_vec);

M1 = moment_vec(2);
M2 = moment_vec(11);

% Find the strain of the most tensioned/compressed fiber due to bending
% The height of the scissor member is estimated to be 0.16 m based on the 
% provided figure. Thus, the distance from neutral axis is 0.08 m. 

strain1=M1/barE_A6N01/barI_A6N01*0.08/10^-6; %(unit ppm)
strain2=M2/barE_A6N01/barI_A6N01*0.08/10^-6; %(unit ppm)

fprintf('The maximum strain at the scissor hinge location are: %d and %d \n', strain1, strain2)

