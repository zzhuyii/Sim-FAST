clear all;
clc;
close all;

%% Initialize the truss
% Define the nodes
H=1; % Height of the truss
L=1; % Length of each span
N=2; % Number of spans

barA=66.8/100/100; % cross section area (m^2)

barE=200*10^9; % Youngs' modulus (Pa)

% We will also solve for the bending stiffness of the frame
I=3124/(100^4);

% When using four segment to represent a frame, 
% three rotational springs are used
% Using Larry Howell's pseudo rigid body model
% We can convert the bending stiffness to spring stiffness
barL=sqrt(H^2+L^2);
kspr=3*barE*I/barL;



%% Define the geometry
% Define the node object
node=Elements_Nodes;
% Create the bar object
bar=Vec_Elements_Bars;
% Create the 3-node rot spring
rot_spr_3N=CD_Elements_RotSprings_3N;

% Create tje truss Assembly
assembly=Assembly_2D_Mechanism;

% This assembly has a node object and a bar object
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_3N=rot_spr_3N;



%% Define the nodal coordinates
% Here we need to set up the nodal coordinates of our truss
node.coordinates_mat=[0*L 0 0; 
                      0*L 0 H;];

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
      (i-1)+L/4  0  H/4;
      (i-1)+L/4  0  3*H/4;
      (i-1)+L/2  0  H/2;
      (i-1)+3*L/4  0  H/4;
      (i-1)+3*L/4  0  3*H/4;

      (i-1)+L/4  0  0;
      (i-1)+L/2  0  0;
      (i-1)+3*L/4  0  0;  % bottom

      (i-1)+L/4  0  H;
      (i-1)+L/2  0  H;
      (i-1)+3*L/4  0  H;  % top

      (i-1)+L  0  0;
      (i-1)+L  0  H;
      ];
end



%% Set up the plotting function for inspection
plots=Plot_2D_Mechanism();
plots.assembly=assembly;

% Range of the plot
plots.displayRange=[-1;5;-1;1;-1;2]; 

% Change the viewing angle
plots.viewAngle1=0;
plots.viewAngle2=0;

% Plot the nodal coordinates for inspection
plots.Plot_Shape_Node_Number()


%% Bar Definition
for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        (i-1)*13+1 (i-1)*13+3;
        (i-1)*13+3 (i-1)*13+5;
        (i-1)*13+5 (i-1)*13+7;
        (i-1)*13+7 (i-1)*13+15;

        (i-1)*13+2 (i-1)*13+4;
        (i-1)*13+4 (i-1)*13+5;
        (i-1)*13+5 (i-1)*13+6;
        (i-1)*13+6 (i-1)*13+14;

        (i-1)*13+1 (i-1)*13+8;
        (i-1)*13+8 (i-1)*13+9;
        (i-1)*13+9 (i-1)*13+10;
        (i-1)*13+10 (i-1)*13+14;

        (i-1)*13+2 (i-1)*13+11;
        (i-1)*13+11 (i-1)*13+12;
        (i-1)*13+12 (i-1)*13+13;
        (i-1)*13+13 (i-1)*13+15;
        ];    

end

% Define the area of the bars
barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);

% Initialize the entire assembly 
assembly.Initialize_Assembly();
% Plot the bar numbering for inspection
plots.Plot_Shape_Bar_Number()




%% Define rotational springs

for i=1:N
    rot_spr_3N.node_ijk_mat=[rot_spr_3N.node_ijk_mat;
        (i-1)*13+1 (i-1)*13+3 (i-1)*13+5;
        (i-1)*13+3 (i-1)*13+5 (i-1)*13+7;
        (i-1)*13+5 (i-1)*13+7 (i-1)*13+13;

        (i-1)*13+2 (i-1)*13+4 (i-1)*13+5;
        (i-1)*13+4 (i-1)*13+5 (i-1)*13+6;
        (i-1)*13+5 (i-1)*13+6 (i-1)*13+14;

        (i-1)*13+1 (i-1)*13+8 (i-1)*13+9;
        (i-1)*13+8 (i-1)*13+9 (i-1)*13+10;
        (i-1)*13+9 (i-1)*13+10 (i-1)*13+14;

        (i-1)*13+2 (i-1)*13+11 (i-1)*13+12;
        (i-1)*13+11 (i-1)*13+12 (i-1)*13+13;
        (i-1)*13+12 (i-1)*13+13 (i-1)*13+15;
        ]; 
end

% Define the rotational spring stiffness
rotNum=size(rot_spr_3N.node_ijk_mat,1);
rot_spr_3N.rot_spr_K_vec=kspr*ones(rotNum,1);


% Initialize the entire assembly again for the new rot-spring
assembly.Initialize_Assembly();

% Plot the rotational spring number
plots.Plot_Shape_Spr_Number





%% Set up the loading solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nodeNum=size(node.coordinates_mat,1);
nr.supp=[(1:nodeNum)', zeros(nodeNum,1), ones(nodeNum,1), zeros(nodeNum,1)];

nr.supp(1,:)=[1,1,1,1];
nr.supp(27,:)=[27,1,1,1];

[F,K]= assembly.Solve_FK(zeros(nodeNum,3));
[Frs3,Krs3]= assembly.rot_spr_3N.Solve_FK(node,zeros(nodeNum,3));
figure
spy(Krs3)

% Set up the load
nr.load=[14,0,0,-100];

% Set up the total loading step
nr.increStep=5;
% Set up the maximum iteration
nr.iterMax=30;
% Set up the tolorence
nr.tol=10^-6;


% Solve for the deformation history
Uhis=nr.Solve();

% Plot the deformed shape
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)),zeros(nodeNum,3));

%% Processing the data
rotSF=rot_spr_3N.theta_stress_free_vec;
rot_spr_3N.Solve_Global_Theta(node,squeeze(Uhis(end,:,:)));
rotC=rot_spr_3N.theta_current_vec;
rotK=rot_spr_3N.rot_spr_K_vec;

RotM=abs(rotK.*(rotC-rotSF));

figure

View1=plots.viewAngle1;
View2=plots.viewAngle2;
view(View1,View2); 

Vsize=plots.displayRange;
axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])

set(gca,'DataAspectRatio',[1 1 1])
set(gcf, 'color', 'white');
set(gcf,'position',[plots.x0,plots.y0,plots.width,plots.height])

barConnect=bar.node_ij_mat;

for j=1:barNum
    node1=assembly.node.coordinates_mat(barConnect(j,1),:);
    node2=assembly.node.coordinates_mat(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color','k');
end

momentFactor=0.0004;

for i=1:N
    p1=node.coordinates_mat((i-1)*13+1,:);
    p2=node.coordinates_mat((i-1)*13+3,:)+[-sqrt(2) 0 sqrt(2)]*RotM((i-1)*12+1)*momentFactor;
    p3=node.coordinates_mat((i-1)*13+5,:)+[-sqrt(2) 0 sqrt(2)]*RotM((i-1)*12+2)*momentFactor;
    p4=node.coordinates_mat((i-1)*13+7,:)+[-sqrt(2) 0 sqrt(2)]*RotM((i-1)*12+3)*momentFactor;
    p5=node.coordinates_mat((i-1)*13+15,:);

    px=[p1(1) p2(1) p3(1) p4(1) p5(1)];
    py=[p1(2) p2(2) p3(2) p4(2) p5(2)];
    pz=[p1(3) p2(3) p3(3) p4(3) p5(3)];
    patch(px,py,pz,'blue','FaceAlpha',.3)


    p1=node.coordinates_mat((i-1)*13+2,:);
    p2=node.coordinates_mat((i-1)*13+4,:)+[sqrt(2) 0 sqrt(2)]*RotM((i-1)*12+4)*momentFactor;
    p3=node.coordinates_mat((i-1)*13+5,:)+[sqrt(2) 0 sqrt(2)]*RotM((i-1)*12+5)*momentFactor;
    p4=node.coordinates_mat((i-1)*13+6,:)+[sqrt(2) 0 sqrt(2)]*RotM((i-1)*12+6)*momentFactor;
    p5=node.coordinates_mat((i-1)*13+14,:);

    px=[p1(1) p2(1) p3(1) p4(1) p5(1)];
    py=[p1(2) p2(2) p3(2) p4(2) p5(2)];
    pz=[p1(3) p2(3) p3(3) p4(3) p5(3)];
    patch(px,py,pz,'blue','FaceAlpha',.3)

end


sigmaY=300*10^6;
MY=sigmaY*I/0.203*2;

Mana=max(RotM);

Fu=5*MY/Mana; %(N)


Lsum=sum(bar.L0_vec);

SelfWeight=Lsum*52.4*9.8;

Fu/SelfWeight
