clear all
close all
clc

tic

%% Define Geometry

sectionNum=3;

layerL=0.06;
holeL=0.004;
layerH=0.002;

t=0.00011; % 0.11mm
E=1.3*10^9;
v=0.2;

kspr=E*t^3/12;
factor=0.14;
barA=t*layerL*factor;

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


% assembly.actBar=actBar;
assembly.cst=cst;
assembly.rotSpr=rotSpr;

%% Define Plotting Functions
plots=Plot_Membrane;
plots.assembly=assembly;
plots.displayRange=[-1; 1; -1; 1; -1; sectionNum/3]*layerL;

plots.Plot_Shape_NodeNumber;


%% Set up the triangle information for pressure
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
    tri_direction=[tri_direction;
        ones(28,1)];
end

triNum=size(tri_ijk,1);

% Organize the trinagle sequence for the rigth direction
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

%% Define Triangle
cst.node_ijk_mat=tri_ijk;
cstNum=size(cst.node_ijk_mat,1);
cst.E_vec=E*ones(cstNum,1);
cst.t_vec=t*ones(cstNum,1);
cst.v_vec=v*ones(cstNum,1);

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

rotSpr.theta1=0.3*pi;
rotSpr.theta2=(2-0.3)*pi;

% Find the bending stiffness
sprNum=size(rotSpr.node_ijkl_mat,1);
rotSpr.rot_spr_K_vec=kspr*ones(sprNum,1);
plots.Plot_Shape_SprNumber;

% Initialize the assembly
assembly.Initialize_Assembly;

% Check force matrix
xcurrent=node.coordinates_mat+node.current_U_mat;
fmat=SolvePressureForce(xcurrent,100,tri_ijk);

%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    1 1 1 1;    
    2 1 1 1;
    3 1 1 1;
    4 1 1 1;
];

step=150;
targetP=20000;


Uhis=[];
for k=1:step

    if k<5
        subStep=50;
        for q=1:subStep
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
        Uhis(k,:,:)=squeeze(nr.Solve());
    end

    xcurrent=node.coordinates_mat+node.current_U_mat;

    nodeNum=size(node.coordinates_mat,1);
    fmat=SolvePressureForce(xcurrent,k/step*targetP,tri_ijk);
    fpressure=[(1:nodeNum)',fmat];

    nr.load=fpressure;
    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=10^-8;

    Uhis(k,:,:)=squeeze(nr.Solve());

end
toc

figure
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))


%% Blocking Force
B_Force=0.5; % Blocking Force
B_dir=[0,1,0]; % Blocking direction

nr.supp=[
    13 1 1 1;    
    15 1 1 1;
    17 1 1 1;
    19 1 1 1;
];
step=200;

Uhis=[];
for k=1:step

    xcurrent=node.coordinates_mat+node.current_U_mat;

    nodeNum=size(node.coordinates_mat,1);
    fmat=SolvePressureForce(xcurrent,targetP,tri_ijk);
    fpressure=[(1:nodeNum)',fmat];

    % First step, put 50 grams on nodes
    fpressure(17+(sectionNum-1)*28,2:4)=fpressure(17+(sectionNum-1)*28,2:4) + B_dir*B_Force*k/step;
    fpressure(19+(sectionNum-1)*28,2:4)=fpressure(19+(sectionNum-1)*28,2:4) + B_dir*B_Force*k/step;

    nr.load=fpressure;
    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=10^-8;

    Uhis(k,:,:)=squeeze(nr.Solve());

end

% Compute Bending Angle
angle_v1 = xcurrent(72,:)-xcurrent(70,:);
angle_v2 = xcurrent(76,:)-xcurrent(70,:);
angle_cross_v = cross(angle_v1,angle_v2);
angle_z = [0,0,1];
cosTheta = dot(angle_cross_v, angle_z) / (norm(angle_cross_v) * norm(angle_z));
theta_rad = acos(cosTheta);
theta_deg_Blocking_initial = rad2deg(theta_rad); % Bending Angle at initial stage

figure
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))

UhisTotal=Uhis;



%% Bending
nr.supp=[
    13 1 1 1;    
    15 1 1 1;
    17 1 1 1;
    19 1 1 1;
];

step=200;
targetF=0.96; % Force on each cable for bending

Uhis=[];
for k=1:step

    xcurrent=node.coordinates_mat+node.current_U_mat;

    nodeNum=size(node.coordinates_mat,1);
    fmat=SolvePressureForce(xcurrent,targetP,tri_ijk);
    fpressure=[(1:nodeNum)',fmat];

    for i=1:sectionNum-1
    
        n1=xcurrent((i)*28+13,:);
        n2=xcurrent((i-1)*28+13,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*28+13,2:4)=fpressure((i)*28+13,2:4)+v1*targetF*k/step;
        fpressure((i-1)*28+13,2:4)=fpressure((i-1)*28+13,2:4)+v2*targetF*k/step;

        n1=xcurrent((i)*28+15,:);
        n2=xcurrent((i-1)*28+15,:);

        v1=n2-n1;
        v2=n1-n2;

        v1=v1/norm(v1);
        v2=v2/norm(v2);

        fpressure((i)*28+15,2:4)=fpressure((i)*28+15,2:4)+v1*targetF*k/step;
        fpressure((i-1)*28+15,2:4)=fpressure((i-1)*28+15,2:4)+v2*targetF*k/step;

    end

    fpressure(17+(sectionNum-1)*28,2:4)=fpressure(17+(sectionNum-1)*28,2:4) + B_dir*B_Force; % 50 grams still applied on the structure
    fpressure(19+(sectionNum-1)*28,2:4)=fpressure(19+(sectionNum-1)*28,2:4) + B_dir*B_Force; % 50 grams still applied on the structure

    nr.load=fpressure;
    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=10^-8;

    Uhis(k,:,:)=squeeze(nr.Solve());

end


% Compute Bending Angle
angle_v1 = xcurrent(72,:)-xcurrent(70,:);
angle_v2 = xcurrent(76,:)-xcurrent(70,:);
angle_cross_v = cross(angle_v1,angle_v2);
angle_z = [0,0,1];
cosTheta = dot(angle_cross_v, angle_z) / (norm(angle_cross_v) * norm(angle_z));
theta_rad = acos(cosTheta);
theta_deg_Blocking = rad2deg(theta_rad); % Bending Angle

figure
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))

figure
UhisTotal=cat(1,UhisTotal,Uhis);
plots.fileName='Bending_Angle_Blocking_Force.gif';
plots.Plot_DeformedHis(UhisTotal(1:10:end,:,:))



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
