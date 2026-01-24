clear all
close all
clc
tic

%% Define Geometry
sectionNum=3;
layerL=0.06;
holeL=0.004;
layerH=0.002;

% film properties
t=0.00011; % 0.11mm for the thickness of film
E=1.3*10^9; % Young's modulus of the film
v=0.2; % Poisson's Ratio of the film

% rotational spring's stiffness
kspr=E*t^3/12; % The rotational springs' stiffness

% target pressure
targetP=10000;

%% Define the Nodal Coordinates
node=Elements_Nodes;
for i=1:sectionNum
    node.coordinates_mat=[node.coordinates_mat
        -holeL/2 -holeL/2 (i-1)*layerH;
        holeL/2 -holeL/2 (i-1)*layerH;
        holeL/2 holeL/2 (i-1)*layerH;
        -holeL/2 holeL/2 (i-1)*layerH;

        -layerL/4 -layerL/4 (i-0.8)*layerH;
        0 -layerL/4 (i-0.8)*layerH;
        layerL/4 -layerL/4 (i-0.8)*layerH;
        layerL/4 0 (i-0.8)*layerH;
        layerL/4 layerL/4 (i-0.85)*layerH;
        0 layerL/4 (i-0.8)*layerH;
        -layerL/4 layerL/4 (i-0.8)*layerH;
        -layerL/4 0 (i-0.8)*layerH;

        -layerL/2 -layerL/2 (i-0.5)*layerH;
        0 -layerL/2 (i-0.5)*layerH;
        layerL/2 -layerL/2 (i-0.5)*layerH;
        layerL/2 0 (i-0.5)*layerH;
        layerL/2 layerL/2 (i-0.5)*layerH;
        0 layerL/2 (i-0.5)*layerH;
        -layerL/2 layerL/2 (i-0.5)*layerH;
        -layerL/2 0 (i-0.5)*layerH;

        -layerL/4 -layerL/4 (i-0.2)*layerH;
        0 -layerL/4 (i-0.2)*layerH;
        layerL/4 -layerL/4 (i-0.2)*layerH;
        layerL/4 0 (i-0.2)*layerH;
        layerL/4 layerL/4 (i-0.2)*layerH;
        0 layerL/4 (i-0.2)*layerH;
        -layerL/4 layerL/4 (i-0.2)*layerH;
        -layerL/4 0 (i-0.2)*layerH;
        ];
end

node.coordinates_mat=[node.coordinates_mat
    -holeL/2 -holeL/2 sectionNum*layerH;
    holeL/2 -holeL/2 sectionNum*layerH;
    holeL/2 holeL/2 sectionNum*layerH;
    -holeL/2 holeL/2 sectionNum*layerH];


%% Define assembly
assembly=Assembly_Membrane;
assembly.node=node;

cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings_4N;

assembly.cst=cst;
assembly.rotSpr=rotSpr;

%% Define Plotting Functions
plots=Plot_Membrane;
plots.assembly=assembly;
plots.displayRange=[-1; 1; -1; 1; -1; sectionNum/3]*layerL;

% plot the nodal location for inspection
plots.Plot_Shape_NodeNumber;


%% Set up the Triangle Elements
tri_ijk=[];
tri_direction=[];

for i=1:sectionNum  
    tri_ijk=[tri_ijk;
        (i-1)*28+1 (i-1)*28+2 (i-1)*28+6;    
        (i-1)*28+2 (i-1)*28+3 (i-1)*28+8;   
        (i-1)*28+3 (i-1)*28+4 (i-1)*28+10;   
        (i-1)*28+4 (i-1)*28+1 (i-1)*28+12;

        (i-1)*28+1 (i-1)*28+5 (i-1)*28+6;    
        (i-1)*28+2 (i-1)*28+6 (i-1)*28+7;   
        (i-1)*28+2 (i-1)*28+7 (i-1)*28+8;   
        (i-1)*28+3 (i-1)*28+8 (i-1)*28+9;
        (i-1)*28+3 (i-1)*28+9 (i-1)*28+10; 
        (i-1)*28+4 (i-1)*28+10 (i-1)*28+11; 
        (i-1)*28+4 (i-1)*28+11 (i-1)*28+12;
        (i-1)*28+1 (i-1)*28+12 (i-1)*28+5; 

        (i-1)*28+5 (i-1)*28+13 (i-1)*28+14;    
        (i-1)*28+5 (i-1)*28+6 (i-1)*28+14;   
        (i-1)*28+6 (i-1)*28+7 (i-1)*28+14;   
        (i-1)*28+7 (i-1)*28+14 (i-1)*28+15;
        (i-1)*28+7 (i-1)*28+15 (i-1)*28+16; 
        (i-1)*28+7 (i-1)*28+8 (i-1)*28+16; 
        (i-1)*28+8 (i-1)*28+9 (i-1)*28+16;
        (i-1)*28+9 (i-1)*28+16 (i-1)*28+17; 

        (i-1)*28+9 (i-1)*28+17 (i-1)*28+18;    
        (i-1)*28+9 (i-1)*28+10 (i-1)*28+18;   
        (i-1)*28+10 (i-1)*28+11 (i-1)*28+18;   
        (i-1)*28+11 (i-1)*28+18 (i-1)*28+19;
        (i-1)*28+11 (i-1)*28+19 (i-1)*28+20; 
        (i-1)*28+11 (i-1)*28+12 (i-1)*28+20; 
        (i-1)*28+12 (i-1)*28+5 (i-1)*28+20;
        (i-1)*28+5 (i-1)*28+20 (i-1)*28+13; 
        ];   

    % The bottom film is facing down
    tri_direction=[tri_direction;
        -ones(28,1)];

    tri_ijk=[tri_ijk;
        (i)*28+1 (i)*28+2 (i-1)*28+6+16;    
        (i)*28+2 (i)*28+3 (i-1)*28+8+16;   
        (i)*28+3 (i)*28+4 (i-1)*28+10+16;   
        (i)*28+4 (i)*28+1 (i-1)*28+12+16;

        (i)*28+1 (i-1)*28+5+16 (i-1)*28+6+16;    
        (i)*28+2 (i-1)*28+6+16 (i-1)*28+7+16;   
        (i)*28+2 (i-1)*28+7+16 (i-1)*28+8+16;   
        (i)*28+3 (i-1)*28+8+16 (i-1)*28+9+16;
        (i)*28+3 (i-1)*28+9+16 (i-1)*28+10+16; 
        (i)*28+4 (i-1)*28+10+16 (i-1)*28+11+16; 
        (i)*28+4 (i-1)*28+11+16 (i-1)*28+12+16;
        (i)*28+1 (i-1)*28+12+16 (i-1)*28+5+16; 

        (i-1)*28+5+16 (i-1)*28+13 (i-1)*28+14;    
        (i-1)*28+5+16 (i-1)*28+6+16 (i-1)*28+14;   
        (i-1)*28+6+16 (i-1)*28+7+16 (i-1)*28+14;   
        (i-1)*28+7+16 (i-1)*28+14 (i-1)*28+15;
        (i-1)*28+7+16 (i-1)*28+15 (i-1)*28+16; 
        (i-1)*28+7+16 (i-1)*28+8+16 (i-1)*28+16; 
        (i-1)*28+8+16 (i-1)*28+9+16 (i-1)*28+16;
        (i-1)*28+9+16 (i-1)*28+16 (i-1)*28+17; 

        (i-1)*28+9+16 (i-1)*28+17 (i-1)*28+18;    
        (i-1)*28+9+16 (i-1)*28+10+16 (i-1)*28+18;   
        (i-1)*28+10+16 (i-1)*28+11+16 (i-1)*28+18;   
        (i-1)*28+11+16 (i-1)*28+18 (i-1)*28+19;
        (i-1)*28+11+16 (i-1)*28+19 (i-1)*28+20; 
        (i-1)*28+11+16 (i-1)*28+12+16 (i-1)*28+20; 
        (i-1)*28+12+16 (i-1)*28+5+16 (i-1)*28+20;
        (i-1)*28+5+16 (i-1)*28+20 (i-1)*28+13; 
        ];   

    % The top film is facing up
    tri_direction=[tri_direction;
        ones(28,1)];
end

% number of triangles
triNum=size(tri_ijk,1);

% Organize the trinagle sequence so that 
% the normal vector points to the out side
for i=1:triNum
    x1=node.coordinates_mat(tri_ijk(i,1),:);
    x2=node.coordinates_mat(tri_ijk(i,2),:);
    x3=node.coordinates_mat(tri_ijk(i,3),:);
    v1=x2-x1;
    v2=x3-x1;
    normal=cross(v1,v2);
    if sign(normal(3))==tri_direction(i)
    else
        temp=tri_ijk(i,2);
        tri_ijk(i,2)=tri_ijk(i,3);
        tri_ijk(i,3)=temp;
    end
end

% Define Triangle Elements
cst.node_ijk_mat=tri_ijk;

% Define Young's modulus, thickness, Poisson's ratio
cst.E_vec=E*ones(triNum,1);
cst.t_vec=t*ones(triNum,1);
cst.v_vec=v*ones(triNum,1);

% Plot the CST elements for inspection
plots.Plot_Shape_CSTNumber;

%% Define Rotational Spring
for i=1:sectionNum
    rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;

        (i-1)*28+1 (i-1)*28+5 (i-1)*28+6 (i-1)*28+14;
        (i-1)*28+2 (i-1)*28+6 (i-1)*28+7 (i-1)*28+14;
        (i-1)*28+2 (i-1)*28+7 (i-1)*28+8 (i-1)*28+16;
        (i-1)*28+3 (i-1)*28+8 (i-1)*28+9 (i-1)*28+16;
        (i-1)*28+3 (i-1)*28+9 (i-1)*28+10 (i-1)*28+18;
        (i-1)*28+4 (i-1)*28+10 (i-1)*28+11 (i-1)*28+18;
        (i-1)*28+4 (i-1)*28+11 (i-1)*28+12 (i-1)*28+20;
        (i-1)*28+1 (i-1)*28+12 (i-1)*28+5 (i-1)*28+20;

        (i)*28+1 (i-1)*28+5+16 (i-1)*28+6+16 (i-1)*28+14;
        (i)*28+2 (i-1)*28+6+16 (i-1)*28+7+16 (i-1)*28+14;
        (i)*28+2 (i-1)*28+7+16 (i-1)*28+8+16 (i-1)*28+16;
        (i)*28+3 (i-1)*28+8+16 (i-1)*28+9+16 (i-1)*28+16;
        (i)*28+3 (i-1)*28+9+16 (i-1)*28+10+16 (i-1)*28+18;
        (i)*28+4 (i-1)*28+10+16 (i-1)*28+11+16 (i-1)*28+18;
        (i)*28+4 (i-1)*28+11+16 (i-1)*28+12+16 (i-1)*28+20;
        (i)*28+1 (i-1)*28+12+16 (i-1)*28+5+16 (i-1)*28+20;

        (i-1)*28+5 (i-1)*28+14 (i-1)*28+6 (i-1)*28+7;
        (i-1)*28+7 (i-1)*28+16 (i-1)*28+8 (i-1)*28+9;
        (i-1)*28+9 (i-1)*28+18 (i-1)*28+10 (i-1)*28+11;
        (i-1)*28+11 (i-1)*28+20 (i-1)*28+12 (i-1)*28+5;

        (i-1)*28+5+16 (i-1)*28+14 (i-1)*28+6+16 (i-1)*28+7+16;
        (i-1)*28+7+16 (i-1)*28+16 (i-1)*28+8+16 (i-1)*28+9+16;
        (i-1)*28+9+16 (i-1)*28+18 (i-1)*28+10+16 (i-1)*28+11+16;
        (i-1)*28+11+16 (i-1)*28+20 (i-1)*28+12+16 (i-1)*28+5+16;

        (i-1)*28+5 (i-1)*28+1 (i-1)*28+6 (i-1)*28+2;
        (i-1)*28+7 (i-1)*28+2 (i-1)*28+8 (i-1)*28+3;
        (i-1)*28+9 (i-1)*28+3 (i-1)*28+10 (i-1)*28+4;
        (i-1)*28+11 (i-1)*28+4 (i-1)*28+12 (i-1)*28+1;

        (i-1)*28+1 (i-1)*28+2 (i-1)*28+6 (i-1)*28+7;
        (i-1)*28+2 (i-1)*28+3 (i-1)*28+8 (i-1)*28+9;
        (i-1)*28+3 (i-1)*28+4 (i-1)*28+10 (i-1)*28+11;
        (i-1)*28+4 (i-1)*28+1 (i-1)*28+12 (i-1)*28+5;

        (i-1)*28+5+16 (i)*28+1 (i-1)*28+6+16 (i)*28+2;
        (i-1)*28+7+16 (i)*28+2 (i-1)*28+8+16 (i)*28+3;
        (i-1)*28+9+16 (i)*28+3 (i-1)*28+10+16 (i)*28+4;
        (i-1)*28+11+16 (i)*28+4 (i-1)*28+12+16 (i)*28+1;

        (i)*28+1 (i)*28+2 (i-1)*28+6+16 (i-1)*28+7+16;
        (i)*28+2 (i)*28+3 (i-1)*28+8+16 (i-1)*28+9+16;
        (i)*28+3 (i)*28+4 (i-1)*28+10+16 (i-1)*28+11+16;
        (i)*28+4 (i)*28+1 (i-1)*28+12+16 (i-1)*28+5+16;
        ];

    if i<sectionNum
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
            (i-1)*28+22 (i-1)*28+29 (i-1)*28+30 (i-1)*28+34;
            (i-1)*28+24 (i-1)*28+30 (i-1)*28+31 (i-1)*28+36;
            (i-1)*28+26 (i-1)*28+31 (i-1)*28+32 (i-1)*28+38;
            (i-1)*28+28 (i-1)*28+32 (i-1)*28+29 (i-1)*28+40;
            ];
    end
end

% Rotational angle where contact will come into effects
rotSpr.theta1=0.1*pi;
rotSpr.theta2=(2-0.1)*pi;

% Define the bending stiffness
sprNum=size(rotSpr.node_ijkl_mat,1);
rotSpr.rot_spr_K_vec=kspr*ones(sprNum,1);

% Plot the rotational springs for inspection
plots.Plot_Shape_SprNumber;


%% Initialize the assembly
assembly.Initialize_Assembly;


%% Set up solver for inflation process
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    1 1 1 1;    
    2 1 1 1;
    3 1 1 1;
    4 1 1 1;
];

step=50;

UhisInflate=[];
for k=1:step
    % To help with convergence, we further make the first few steps 
    % to have smaller step lengths 
    if k<2
        subStep=10;
        for q=1:subStep
            % Update the pressure forces for nonlinearity
            xcurrent=node.coordinates_mat+node.current_U_mat;
            nodeNum=size(node.coordinates_mat,1);
            fmat=SolvePressureForce(xcurrent,(k-1+q/subStep)/step*targetP,tri_ijk);
            fpressure=[(1:nodeNum)',fmat];

            nr.load=fpressure;
            nr.increStep=1;
            nr.iterMax=20;
            nr.tol=10^-8;

            Utemp=nr.Solve();
        end
        UhisInflate(k,:,:)=squeeze(nr.Solve());
    end

    % Update the pressure forces for nonlinearity
    xcurrent=node.coordinates_mat+node.current_U_mat;
    nodeNum=size(node.coordinates_mat,1);
    fmat=SolvePressureForce(xcurrent,k/step*targetP,tri_ijk);
    fpressure=[(1:nodeNum)',fmat];

    nr.load=fpressure;    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=10^-6;

    UhisInflate(k,:,:)=squeeze(nr.Solve());

end

toc
% Plot the inflated configuration
figure;
plots.Plot_DeformedShape(squeeze(UhisInflate(end,:,:)))
plots.fileName='Three_Unit_Compare.gif';
plots.Plot_DeformedHis(UhisInflate);




%% Solve the applied force due to pressure
function Fmat=SolvePressureForce(NodeCord,P,tri_ijk)
    nodeNum=size(NodeCord,1);
    Fmat=zeros(nodeNum,3);
    triNum=size(tri_ijk,1);

    for i=1:triNum
        x1=NodeCord(tri_ijk(i,1),:);
        x2=NodeCord(tri_ijk(i,2),:);
        x3=NodeCord(tri_ijk(i,3),:);

        v1=x2-x1;
        v2=x3-x1;

        normal=cross(v1,v2);
        A=norm(normal)/2;

        f=A*P/3;
        normal=normal/norm(normal);

        Fmat(tri_ijk(i,1),:)=Fmat(tri_ijk(i,1),:)+f*normal;
        Fmat(tri_ijk(i,2),:)=Fmat(tri_ijk(i,2),:)+f*normal;
        Fmat(tri_ijk(i,3),:)=Fmat(tri_ijk(i,3),:)+f*normal;
    end   
end
