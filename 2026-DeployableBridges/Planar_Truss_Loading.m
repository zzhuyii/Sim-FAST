clear all; 
close all; 
clc;

%% Initialize the truss
% Define the nodes
H=1; % Height of the truss (m)
L=1; % Length of each span (m)
N=4;
barA=66.8/100/100; % cross section area (m^2)
barE=200*10^9; % Youngs' modulus (Pa)



%% Define the assembly

node=Elements_Nodes;
bar=CD_Elements_Bars;

assembly=Assembly_2D_Truss;
assembly.node=node;
assembly.bar=bar;



%% Define the nodal coordinates
% Here we need to set up the nodal coordinates of our truss
node.coordinates_mat=[0*L 0 0; 
                      0*L 0 H];


for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
                      (i-0.5)*L 0 0; 
                      (i-0.5)*L 0 0.5*H; 
                      (i-0.5)*L 0 H
                      (i)*L 0 0;  
                      (i)*L 0 H];
end



%% Set up the plotting function for inspection
plots=Plot_2D_Truss();
plots.assembly=assembly;

plots.displayRange=[-1;5;-1;1;-1;2]; 
plots.viewAngle1=0;
plots.viewAngle2=0;

plots.Plot_Shape_Node_Number()


%% Define bars
bar.node_ij_mat=[1 2];
for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        (i-1)*5+1,(i-1)*5+4;
        (i-1)*5+4,(i-1)*5+7;
        (i-1)*5+2,(i-1)*5+4;
        (i-1)*5+4,(i-1)*5+6;
        (i-1)*5+3,(i-1)*5+4;
        (i-1)*5+4,(i-1)*5+5;
        (i-1)*5+1,(i-1)*5+3;
        (i-1)*5+3,(i-1)*5+6;
        (i-1)*5+2,(i-1)*5+5;
        (i-1)*5+5,(i-1)*5+7;
        (i-1)*5+6,(i-1)*5+7;]; 

end

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);


assembly.Initialize_Assembly();
plots.Plot_Shape_Bar_Number()


%% Set up the loading solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nodeNum=size(node.coordinates_mat,1);
nr.supp=[(1:nodeNum)', zeros(nodeNum,1), ones(nodeNum,1), zeros(nodeNum,1)];

nr.supp(1,:)=[1,1,1,1];
nr.supp(N*5+1,:)=[N*5+1,1,1,1];

% Set up the load
nr.load=[N*5/2+1,0,0,-100];

% Set up the total loading step
nr.increStep=10;
% Set up the maximum iteration
nr.iterMax=30;

% Solve for the deformation history
Uhis=nr.Solve();

% Plot the deformed shape
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)));
truss_strain=bar.Solve_Strain(node,squeeze(Uhis(end,:,:)));
internal_force=(truss_strain).*(bar.E_vec).*(bar.A_vec);

plots.Plot_Bar_Force(internal_force);

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

forceFactor=1/5000;

p1=node.coordinates_mat(N*5+1,:);
p2=node.coordinates_mat(N*5+1,:)+[1 0 0]*internal_force(end)*forceFactor;
p3=node.coordinates_mat(N*5+2,:)+[1 0 0]*internal_force(end)*forceFactor;
p4=node.coordinates_mat(N*5+2,:);

px=[p1(1) p2(1) p3(1) p4(1) ];
py=[p1(2) p2(2) p3(2) p4(2) ];
pz=[p1(3) p2(3) p3(3) p4(3) ];
patch(px,py,pz,'blue','FaceAlpha',.3)

for i=1:N
    p1=node.coordinates_mat((i-1)*5+1,:);
    p2=node.coordinates_mat((i-1)*5+1,:)+[1 0 0]*internal_force((i-1)*11+1)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+2,:)+[1 0 0]*internal_force((i-1)*11+1)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+2,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+1,:);
    p2=node.coordinates_mat((i-1)*5+1,:)+[-sqrt(2) 0 sqrt(2)]*internal_force((i-1)*11+2)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+4,:)+[-sqrt(2) 0 sqrt(2)]*internal_force((i-1)*11+2)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+4,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+7,:);
    p2=node.coordinates_mat((i-1)*5+7,:)+[-sqrt(2) 0 sqrt(2)]*internal_force((i-1)*11+3)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+4,:)+[-sqrt(2) 0 sqrt(2)]*internal_force((i-1)*11+3)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+4,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+2,:);
    p2=node.coordinates_mat((i-1)*5+2,:)+[sqrt(2) 0 sqrt(2)]*internal_force((i-1)*11+4)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+4,:)+[sqrt(2) 0 sqrt(2)]*internal_force((i-1)*11+4)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+4,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)
    
    p1=node.coordinates_mat((i-1)*5+6,:);
    p2=node.coordinates_mat((i-1)*5+6,:)+[sqrt(2) 0 sqrt(2)]*internal_force((i-1)*11+5)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+4,:)+[sqrt(2) 0 sqrt(2)]*internal_force((i-1)*11+5)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+4,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+3,:);
    p2=node.coordinates_mat((i-1)*5+3,:)+[1 0 0]*internal_force((i-1)*11+6)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+4,:)+[1 0 0]*internal_force((i-1)*11+6)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+4,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+5,:);
    p2=node.coordinates_mat((i-1)*5+5,:)+[1 0 0]*internal_force((i-1)*11+7)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+4,:)+[1 0 0]*internal_force((i-1)*11+7)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+4,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+1,:);
    p2=node.coordinates_mat((i-1)*5+1,:)+[0 0 1]*internal_force((i-1)*11+8)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+3,:)+[0 0 1]*internal_force((i-1)*11+8)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+3,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+3,:);
    p2=node.coordinates_mat((i-1)*5+3,:)+[0 0 1]*internal_force((i-1)*11+9)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+6,:)+[0 0 1]*internal_force((i-1)*11+9)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+6,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+2,:);
    p2=node.coordinates_mat((i-1)*5+2,:)+[0 0 1]*internal_force((i-1)*11+10)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+5,:)+[0 0 1]*internal_force((i-1)*11+10)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+5,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)

    p1=node.coordinates_mat((i-1)*5+7,:);
    p2=node.coordinates_mat((i-1)*5+7,:)+[0 0 1]*internal_force((i-1)*11+11)*forceFactor;
    p3=node.coordinates_mat((i-1)*5+5,:)+[0 0 1]*internal_force((i-1)*11+11)*forceFactor;
    p4=node.coordinates_mat((i-1)*5+5,:);

    px=[p1(1) p2(1) p3(1) p4(1) ];
    py=[p1(2) p2(2) p3(2) p4(2) ];
    pz=[p1(3) p2(3) p3(3) p4(3) ];
    patch(px,py,pz,'blue','FaceAlpha',.3)
end


sigmaY=300*10^6;
barFY=sigmaY*barA;

barFana=max(internal_force);

Fu=10*100*barFY/barFana; %(N)


Lsum=sum(bar.L0_vec);
SelfWeight=Lsum*52.4*9.8;

Fu/SelfWeight