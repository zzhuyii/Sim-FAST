%% Initialize the solver
clear all;
clc;
close all;
tic


%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% Definition of Geometry
L1=300*10^-6; % Panel side overhang
L2=100*10^-6; % Beam width
L3=600*10^-6; % Panel middle width
L4=400*10^-6; % Panel length
L5=440*10^-6; % Beam Length

BeamSegment=5;

% actuator length. Because the end of beam is covered by Al, they will not
% fold as the PZT is just on the neutal axis. We need to take out these
% contribution so the active folding beam length is 340 um.
La=340*10^-6; 

% Damping for the system
DampingFactorA=90;
DampingFactorB=1*10^-7;

% Sweep Time Set Up
dt=4*10^-7;
totalSweepTime=0.02;

% Sweep frequency set up
f1=200;
f2=700;

% Folding Informaiton
InitialFoldAngleLeft=45;
InitialFoldAngleRight=45;
MirrorBentAngle=1;

Vdc=5; % dc offset voltage
Vpp=10; % peak to peak voltage
d31=-100*10^-12; % pm/V


%% solve for the area of bars 
% Layer information of the PZT beam
tPtTop=0.1*10^-6;
tPZT=1*10^-6;
tPtBTM=0.15*10^-6;
tSiO2=1*10^-6;

tAl=1*10^-6;

ESiO2=66*10^9;
EPt=168*10^9;
EPZT=70*10^9;
EAl=70*10^9;

Eeq=100*10^9;
% Assume an equivalent Young's modulus of 100 GPa

tBeamTotal=tSiO2+tPtTop+tPZT+tPtBTM;
EAbeam=L2*(tSiO2*ESiO2+tPtTop*EPt+tPZT*EPZT+tPtBTM*EPt);
tave=EAbeam/L2/Eeq;

% Assume material is homogenized with 100 GPa
AbeamLong=EAbeam/Eeq/2;
AbeamHor=tave*L5/6;
AbeamDiag=7/4/sqrt(2)*tave*L2;

% Calculation of panel area
A=(L3+2*L1+2*L2)*L4;
tpanel=25*10^-6;
Lsum=2*(L3+2*L1+2*L2)+4*L4+2*sqrt(L4^2+L3^3)+...
    2*sqrt(L4^2+L2^3)+2*sqrt(L4^2+L1^3);
Apanel=2*A*tpanel/Lsum;


%% Solve for the rotational springs stiffness for bending 
AreaMomentBeam=(tSiO2*L2*(tSiO2/2)*ESiO2)+...
    (tPtBTM*L2*(tPtBTM/2+tSiO2)*EPt)+...
    (tPZT*L2*(tPZT/2+tPtBTM+tSiO2)*EPZT)+...
    (tPtTop*L2*(tPtTop/2+tPZT+tPtBTM+tSiO2)*EPt);
c=AreaMomentBeam/EAbeam;

EIcomp=L2*(ESiO2*(1/12*tSiO2^3+tSiO2*(c-tSiO2/2)^2)+...
    EPt*(1/12*tPtBTM^3+tPtBTM*(c-tPtBTM/2-tSiO2)^2)+...
    EPZT*(1/12*tPZT^3+tPZT*(c-tPZT/2-tPtBTM-tSiO2)^2)+...
    EPt*(1/12*tPtTop^3+tPtTop*(c-tPtTop/2-tPZT-tPtBTM-tSiO2)^2));

SprBeam=(BeamSegment+1)*EIcomp/L5;
SprBeamDiag=2*sqrt(2)*EIcomp/L2;
SprPanel=100*SprBeam;


%% Convert voltage to folding angle

DcMomment=d31*Vdc/tPZT*EPZT*L2*tPZT*(c-tPZT/2-tPtBTM-tSiO2);
DcAngle=DcMomment/EIcomp*La*180/pi; % degree

PPMomment=d31*Vpp/tPZT*EPZT*L2*tPZT*(c-tPZT/2-tPtBTM-tSiO2);
PPAngle=PPMomment/EIcomp*La*180/pi;
ExciatationAngle=PPAngle/2; % degree

bar=Vec_Elements_Bars;
rotSpr=Vec_Elements_RotSprings_4N;
node=Elements_Nodes;

%% Add node of the micro-mirror
% Here we define the nodal coordinates

node.coordinates_mat(1,:)=[-L5/3,L3/2+L2,0];
node.coordinates_mat(BeamSegment*2+2+2,:)=[-L5/3,-L3/2-L2,0];

for i=1:BeamSegment+1
    node.coordinates_mat(1+i,:)=[L5*(i-1)/BeamSegment,L3/2+L2,0];
    node.coordinates_mat((BeamSegment+1)+1+i,:)=[L5*(i-1)/BeamSegment,L3/2,0];

    node.coordinates_mat((BeamSegment+1)*2+2+i,:)=[L5*(i-1)/BeamSegment,-L3/2-L2,0];
    node.coordinates_mat((BeamSegment+1)*3+2+i,:)=[L5*(i-1)/BeamSegment,-L3/2,0];
end

node.coordinates_mat((BeamSegment+1)*4+2+1,:)=[L5,L3/2+L2+L1,0];
node.coordinates_mat((BeamSegment+1)*4+2+2,:)=[L5,-L3/2-L2-L1,0];

node.coordinates_mat((BeamSegment+1)*4+2+3,:)=[L5+L4/2,0,0];

node.coordinates_mat((BeamSegment+1)*4+2+4,:)=[L5+L4,L3/2+L2+L1,0];
node.coordinates_mat((BeamSegment+1)*4+2+5,:)=[L5+L4,L3/2,0];
node.coordinates_mat((BeamSegment+1)*4+2+6,:)=[L5+L4,-L3/2,0];
node.coordinates_mat((BeamSegment+1)*4+2+7,:)=[L5+L4,-L3/2-L2-L1,0];



%% Define longitudinal bar elements

bar.node_ij_mat=[1,2];
bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+2,2*(BeamSegment+1)+3];

bar.node_ij_mat=[bar.node_ij_mat; 1,2+(BeamSegment+1)];
bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+2,3*(BeamSegment+1)+2+1];

for i=1:BeamSegment

    bar.node_ij_mat=[bar.node_ij_mat; i+1,i+2];    
    bar.node_ij_mat=[bar.node_ij_mat; BeamSegment+1+1+i,BeamSegment+1+1+i+1];    
    bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+2+i,2*(BeamSegment+1)+3+i];    
    bar.node_ij_mat=[bar.node_ij_mat; 3*(BeamSegment+1)+2+i,3*(BeamSegment+1)+3+i];

end

bar.E_vec=[bar.E_vec;Eeq*ones((BeamSegment+1)*4,1)];
bar.A_vec=[bar.A_vec;AbeamLong*ones((BeamSegment+1)*4,1)];
bar.L0_vec=[bar.L0_vec;L5/BeamSegment*ones((BeamSegment+1)*4,1)];


%% Define horizontal bar elements
for i=1:BeamSegment+1
    bar.node_ij_mat=[bar.node_ij_mat; 1+i,(BeamSegment+1)+1+i];    
    bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+2+i,3*(BeamSegment+1)+2+i];
end

bar.E_vec=[bar.E_vec;Eeq*ones(BeamSegment*2+2,1)];
bar.A_vec=[bar.A_vec;AbeamHor*ones(BeamSegment*2+2,1)];
bar.L0_vec=[bar.L0_vec;L2*ones(BeamSegment*2+2,1)];


%% Define horizontal rotational springs

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat; 1,2,...
    (BeamSegment+1)+1+1,(BeamSegment+1)+2+1];
rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;2*(BeamSegment+1)+1+1,...
    2*(BeamSegment+1)+2+1,3*(BeamSegment+1)+2+1,3*(BeamSegment+1)+3+1];

for i=2:BeamSegment

    if mod(i,2)==1
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat; (BeamSegment+1)+i,1+i,...
            (BeamSegment+1)+1+i,(BeamSegment+1)+2+i];
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;3*(BeamSegment+1)+1+i,...
            2*(BeamSegment+1)+2+i,3*(BeamSegment+1)+2+i,3*(BeamSegment+1)+3+i];
    else
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat; i,1+i,...
            (BeamSegment+1)+1+i,2+i];
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;2*(BeamSegment+1)+1+i,...
            2*(BeamSegment+1)+2+i,3*(BeamSegment+1)+2+i,2*(BeamSegment+1)+3+i];

    end

end

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat; ...
    (BeamSegment+1)+1+BeamSegment,1+BeamSegment+1,...
    (BeamSegment+1)+1+BeamSegment+1,4*(BeamSegment+1)+2+5];
rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    3*(BeamSegment+1)+1+BeamSegment+1,2*(BeamSegment+1)+2+BeamSegment+1,...
    3*(BeamSegment+1)+2+BeamSegment+1,4*(BeamSegment+1)+2+6];

rotSpr.rot_spr_K_vec=[rotSpr.rot_spr_K_vec;SprBeam*ones((BeamSegment+1)*2,1)];

baseSprNum=[1,2,2*BeamSegment+1,2*BeamSegment+2];
rotSpr.rot_spr_K_vec(baseSprNum)=rotSpr.rot_spr_K_vec(baseSprNum)*2;

rotSpr.theta_stress_free_vec=[rotSpr.theta_stress_free_vec;pi*ones((BeamSegment+1)*2,1)];


%% Define diagonal bars

for i=1:BeamSegment
    if mod(i,2)==1
        bar.node_ij_mat=[bar.node_ij_mat; i+1,(BeamSegment+1)+2+i];
        bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+2+i,...
            3*(BeamSegment+1)+3+i];
    else
        bar.node_ij_mat=[bar.node_ij_mat; i+2,(BeamSegment+1)+1+i];
        bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+3+i,...
            3*(BeamSegment+1)+2+i];
    end

end
bar.E_vec=[bar.E_vec;Eeq*ones(BeamSegment*2,1)];
bar.A_vec=[bar.A_vec;AbeamDiag*ones(BeamSegment*2,1)];
bar.L0_vec=[bar.L0_vec;sqrt(L2*L2+L5*L5/9)*ones(BeamSegment*2,1)];


%% Define diagonal rotational springs
for i=1:BeamSegment
    if mod(i,2)==1
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
            (BeamSegment+1)+1+i,i+1,...
            (BeamSegment+1)+2+i,i+2];

        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
            3*(BeamSegment+1)+2+i,2*(BeamSegment+1)+2+i,...
            3*(BeamSegment+1)+3+i,2*(BeamSegment+1)+3+i];
    else
        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;i+1,...
            (BeamSegment+1)+1+i,i+2,(BeamSegment+1)+2+i];

        rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
            2*(BeamSegment+1)+2+i,3*(BeamSegment+1)+2+i,...
            2*(BeamSegment+1)+3+i,3*(BeamSegment+1)+3+i];
    end

end

rotSpr.rot_spr_K_vec=[rotSpr.rot_spr_K_vec;SprBeamDiag*ones(BeamSegment*2,1)];
rotSpr.theta_stress_free_vec=[rotSpr.theta_stress_free_vec;pi*ones(BeamSegment*2,1)];


%% Define bar for the panels (perimeter)
bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+1,4*(BeamSegment+1)+2];
bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2+5,4*(BeamSegment+1)+2+6];

bar.E_vec=[bar.E_vec;Eeq*ones(2,1)];
bar.A_vec=[bar.A_vec;Apanel*ones(2,1)];
bar.L0_vec=[bar.L0_vec;L3*ones(2,1)];

bar.node_ij_mat=[bar.node_ij_mat; (BeamSegment+1)+1,4*(BeamSegment+1)+2+1];
bar.node_ij_mat=[bar.node_ij_mat; 3*(BeamSegment+1)+2,4*(BeamSegment+1)+2+2];

bar.E_vec=[bar.E_vec;Eeq*ones(2,1)];
bar.A_vec=[bar.A_vec;Apanel*ones(2,1)];
bar.L0_vec=[bar.L0_vec;L1*ones(2,1)];

bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2+1,4*(BeamSegment+1)+2+4];
bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2+2,4*(BeamSegment+1)+2+7];

bar.E_vec=[bar.E_vec;Eeq*ones(2,1)];
bar.A_vec=[bar.A_vec;Apanel*ones(2,1)];
bar.L0_vec=[bar.L0_vec;L4*ones(2,1)];

bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2+4,4*(BeamSegment+1)+2+5];
bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2+6,4*(BeamSegment+1)+2+7];

bar.E_vec=[bar.E_vec;Eeq*ones(2,1)];
bar.A_vec=[bar.A_vec;Apanel*ones(2,1)];
bar.L0_vec=[bar.L0_vec;(L1+L2)*ones(2,1)];


%% Define bar for inside the panel
bar.node_ij_mat=[bar.node_ij_mat; BeamSegment+2,4*(BeamSegment+1)+2+4];
bar.node_ij_mat=[bar.node_ij_mat; 3*(BeamSegment+1)+2,4*(BeamSegment+1)+2+7];

bar.E_vec=[bar.E_vec;Eeq*ones(2,1)];
bar.A_vec=[bar.A_vec;Apanel*ones(2,1)];
bar.L0_vec=[bar.L0_vec;sqrt(L1*L1+L4*L4)*ones(2,1)];

bar.node_ij_mat=[bar.node_ij_mat; BeamSegment+2,4*(BeamSegment+1)+2+5];
bar.node_ij_mat=[bar.node_ij_mat; 3*(BeamSegment+1)+2,4*(BeamSegment+1)+2+6];

bar.E_vec=[bar.E_vec;Eeq*ones(2,1)];
bar.A_vec=[bar.A_vec;Apanel*ones(2,1)];
bar.L0_vec=[bar.L0_vec;sqrt(L2*L2+L4*L4)*ones(2,1)];

bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+1,4*(BeamSegment+1)+2+5];
bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2,4*(BeamSegment+1)+2+6];

bar.E_vec=[bar.E_vec;Eeq*ones(2,1)];
bar.A_vec=[bar.A_vec;Apanel*ones(2,1)];
bar.L0_vec=[bar.L0_vec;L4*ones(2,1)];

bar.node_ij_mat=[bar.node_ij_mat; 2*(BeamSegment+1)+1,4*(BeamSegment+1)+2+3];
bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2,4*(BeamSegment+1)+2+3];
bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2+3,4*(BeamSegment+1)+2+5];
bar.node_ij_mat=[bar.node_ij_mat; 4*(BeamSegment+1)+2+3,4*(BeamSegment+1)+2+6];

bar.E_vec=[bar.E_vec;Eeq*ones(4,1)];
bar.A_vec=[bar.A_vec;Apanel*ones(4,1)];
bar.L0_vec=[bar.L0_vec;sqrt(L4*L4/4+L3*L3/4)*ones(4,1)];

%% Define diagonal rotational springs for panel
rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    4*(BeamSegment+1)+2+1,(BeamSegment+1)+1,...
    4*(BeamSegment+1)+2+4,4*(BeamSegment+1)+2+5];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    4*(BeamSegment+1)+2+2,3*(BeamSegment+1)+2,...
    4*(BeamSegment+1)+2+7,4*(BeamSegment+1)+2+6];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    4*(BeamSegment+1)+2+4,(BeamSegment+1)+1,...
    4*(BeamSegment+1)+2+5,2*(BeamSegment+1)+1];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    4*(BeamSegment+1)+2,4*(BeamSegment+1)+2+6,...
    3*(BeamSegment+1)+2,4*(BeamSegment+1)+2+7];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    (BeamSegment+1)+1,2*(BeamSegment+1)+1,...
    4*(BeamSegment+1)+2+5,4*(BeamSegment+1)+2+3];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    3*(BeamSegment+1)+2,4*(BeamSegment+1)+2,...
    4*(BeamSegment+1)+2+6,4*(BeamSegment+1)+2+3];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    2*(BeamSegment+1)+1,4*(BeamSegment+1)+2+3,...
    4*(BeamSegment+1)+2,4*(BeamSegment+1)+2+6];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    4*(BeamSegment+1)+2,4*(BeamSegment+1)+2+3,...
    4*(BeamSegment+1)+2+6,4*(BeamSegment+1)+2+5];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    4*(BeamSegment+1)+2+6,4*(BeamSegment+1)+2+3,...
    4*(BeamSegment+1)+2+5,2*(BeamSegment+1)+1];

rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;...
    4*(BeamSegment+1)+2+5,4*(BeamSegment+1)+2+3,...
    2*(BeamSegment+1)+1,4*(BeamSegment+1)+2];


rotSpr.rot_spr_K_vec=[rotSpr.rot_spr_K_vec;SprPanel*ones(10,1)];
rotSpr.theta_stress_free_vec=[rotSpr.theta_stress_free_vec;pi*ones(10,1)];


%% Initialize and Plot for investigation

assembly=Assembly_MEMS();
assembly.bar=bar;
assembly.node=node;
assembly.rot_spr_4N=rotSpr;

assembly.InitializeAssembly();

plots=Plot_MEMS();
plots.displayRange=2*10^-3;
plots.displayRangeRatio=0.5;
plots.assembly=assembly;

% plots.Plot_Shape_NodeNumber()
% plots.Plot_Shape_BarNumber()
% plots.Plot_Shape_SprNumber()


%% Assign Mass Properties
% density of materials for the beam structure
rhoPZT=8000;
rhoSiO2=2600;
rhoPt=21000;
rhoAl=2700;

% density of materials for panels
rhoSi=2300;
tSi=25*10^-6;

% Solve the mass
node.mass_vec=zeros(4*(BeamSegment+1)+2+7,1);

% Mass node for beams
totalBeamMass=L5*L2*(tSiO2*rhoSiO2+tPtTop*rhoPt...
    +tPtTop*rhoPt+tPZT*rhoPZT);
node.mass_vec(1:4*(BeamSegment+1)+2)=totalBeamMass/((BeamSegment+1)*2);


% Mass node for panels
totalPanelMass=(L1*2+L2*2+L3)*L4*(rhoSi*tSi+rhoAl*tAl);
node.mass_vec([(BeamSegment+1)+1,2*(BeamSegment+1)+1,3*(BeamSegment+1)+4,...
    4*(BeamSegment+1)+2,4*(BeamSegment+1)+2+1,4*(BeamSegment+1)+2+2. ...
    4*(BeamSegment+1)+2+3,4*(BeamSegment+1)+2+4,4*(BeamSegment+1)+2+5,...
    4*(BeamSegment+1)+2+6,4*(BeamSegment+1)+2+7])=...
    node.mass_vec([(BeamSegment+1)+1,2*(BeamSegment+1)+1,3*(BeamSegment+1)+4,...
    4*(BeamSegment+1)+2,4*(BeamSegment+1)+2+1,4*(BeamSegment+1)+2+2. ...
    4*(BeamSegment+1)+2+3,4*(BeamSegment+1)+2+4,4*(BeamSegment+1)+2+5,...
    4*(BeamSegment+1)+2+6,4*(BeamSegment+1)+2+7])...
    +totalPanelMass/30; % node at edge
node.mass_vec(4*(BeamSegment+1)+2+3)=totalPanelMass*2/3; % node at center

node.current_ext_force_mat=zeros(4*(BeamSegment+1)+2+7,3);
node.current_ext_force_mat(:,3)=-node.mass_vec*9.8/sqrt(2);
node.current_ext_force_mat(:,1)=node.mass_vec*9.8/sqrt(2);

%% Frequency analysis for the mode shape
% Obtain the mass matrix and stiffness matrix

[Fvec,Kmat]=assembly.Solve_FK(assembly.node.current_U_mat);
Mmat=node.FindMassMat;

supp=[1,1,1,1;
      2,1,1,1;
      BeamSegment+3,1,1,1;
      2+2*(BeamSegment+1),1,1,1;
      3+2*(BeamSegment+1),1,1,1;
      3+3*(BeamSegment+1),1,1,1];

[Kadj,Fadj]=Mod_K_For_Supp(Kmat,supp,Fvec);
[Madj]=Mod_M_For_Supp(Mmat,supp);

[Umode,frequencySquared]=eig(Kadj,Madj);
frequencySquared=diag(frequencySquared);
[freq,index]=sort(abs(frequencySquared));

freq1=sqrt(freq(1))/pi/2;
freq2=sqrt(freq(2))/pi/2;
freq3=sqrt(freq(3))/pi/2;

Umode1=Umode(:,index(1));
Usize=size(Umode1,1);

Umode1=reshape(Umode1,3,Usize/3)'/10000;

Umode2=Umode(:,index(2));
Umode2=reshape(Umode2,3,Usize/3)'/10000;

Umode3=Umode(:,index(3));
Umode3=reshape(Umode3,3,Usize/3)'/10000;

% plots.Plot_DeformedShape(zeros(size(Umode1)),Umode1)
% plots.Plot_DeformedShape(zeros(size(Umode1)),Umode2)
% plots.Plot_DeformedShape(zeros(size(Umode1)),Umode3)


%% Statically Fold to an Angle
foldAngleRateLeft=InitialFoldAngleLeft/180/(BeamSegment-1)+1;
foldAngleRateRight=InitialFoldAngleRight/180/(BeamSegment-1)...
    +DcAngle/180/(BeamSegment-1)+1;

sf=Solver_NR_Folding_4N();
sf.assembly=assembly;
sf.supp=supp;

activeCreaseLeft=2*(2:(BeamSegment))-1;
activeCreaseRight=2*(2:(BeamSegment));
mirrorFold=[2*(BeamSegment+1)+2*(BeamSegment)+5,...
            2*(BeamSegment+1)+2*(BeamSegment)+6];

sf.increStep=20;
sf.tol=1*10^-8;

sf.targetRot=rotSpr.theta_current_vec;
sf.targetRot(activeCreaseLeft)=foldAngleRateLeft*pi;
sf.targetRot(activeCreaseRight)=foldAngleRateRight*pi;

sf.targetRot(mirrorFold)=(MirrorBentAngle/180+1)*pi;

Uhis=sf.Solve();
% plots.Plot_DeformedShape(zeros(size(Umode1)),squeeze(Uhis(end,:,:)))

%% Solve Frequency after folding 

[Fvec,Kmat]=assembly.Solve_FK(assembly.node.current_U_mat);
Mmat=node.FindMassMat;


[Kadj,Fadj]=Mod_K_For_Supp(Kmat,supp,Fvec);
[Madj]=Mod_M_For_Supp(Mmat,supp);

[Umode,frequencySquared]=eig(Kadj,Madj);
frequencySquared=diag(frequencySquared);
[freq,index]=sort(abs(frequencySquared));


freq1Fold=sqrt(freq(1))/pi/2;
freq2Fold=sqrt(freq(2))/pi/2;
freq3Fold=sqrt(freq(3))/pi/2;

Umode1Fold=Umode(:,index(1));
Usize=size(Umode1Fold,1);
Umode1Fold=reshape(Umode1Fold,3,Usize/3)'/10000;

Umode2Fold=Umode(:,index(2));
Umode2Fold=reshape(Umode2Fold,3,Usize/3)'/10000;

Umode3Fold=Umode(:,index(3));
Umode3Fold=reshape(Umode3Fold,3,Usize/3)'/10000;

% plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)),...
%     squeeze(Uhis(end,:,:))+Umode1Fold)
% plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)),...
%     squeeze(Uhis(end,:,:))+Umode2Fold)
% plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)),...
%     squeeze(Uhis(end,:,:))+Umode3Fold)

%% Solve for transient dynamic responses
caa=Solver_CAA_Dynamics();
caa.assembly=assembly;

caa.supp=supp;

caa.alpha=DampingFactorA;
caa.beta=DampingFactorB;

omega=freq1*2*pi;
dampingRatio1=caa.alpha/2/omega+caa.beta*omega/2;
fprintf('damping ratio 1 is %d \n', dampingRatio1);
omega=freq3*2*pi;
dampingRatio2=caa.alpha/2/omega+caa.beta*omega/2;
fprintf('damping ratio 2 is %d \n', dampingRatio2);

caa.dt=dt;
step=round(totalSweepTime/caa.dt);

% Sine wave input
time=(1:step)*caa.dt;
sprNum=size(rotSpr.rot_spr_K_vec);
sprNum=sprNum(1);

caa.rotSprTargetAngle=pi*ones(step,sprNum);

% Generate sweep function
y = chirp(time,f1,totalSweepTime,f2,'linear', 90);

% Make the mirror bend a little
caa.rotSprTargetAngle(:,2*(BeamSegment+1)+2*(BeamSegment)+5)=...
    (MirrorBentAngle/180+1)*pi;
caa.rotSprTargetAngle(:,2*(BeamSegment+1)+2*(BeamSegment)+6)=...
    (MirrorBentAngle/180+1)*pi;


% One leg is not actuated in the experiment
for i=1:length(activeCreaseLeft)
    caa.rotSprTargetAngle(:,activeCreaseLeft(i))=...
        (foldAngleRateLeft*pi)';
end

% Actuate one leg to generate twisting motion
for i=1:length(activeCreaseRight)
    caa.rotSprTargetAngle(:,activeCreaseRight(i))=...
        (foldAngleRateRight*pi+ExciatationAngle/(BeamSegment-1)/180*pi*y)';
end

% Applied gravity force
caa.Fext=zeros(step, 4*(BeamSegment+1)+2+7,3);
for i=1:step
    caa.Fext(i,:,:)=node.current_ext_force_mat;
end

Uhis=caa.Solve();
CorpFactor=10;
UhisCorp=Uhis(1:CorpFactor:end,:,:);

%% Plot the twisting angle
angleTwist=zeros(step/CorpFactor,1);
vector=zeros(step/CorpFactor,3);

for i=1:step/CorpFactor
    node1=squeeze(UhisCorp(i,4*(BeamSegment+1)+2+7,:))+node.coordinates_mat(4*(BeamSegment+1)+2+7,:)';
    node2=squeeze(UhisCorp(i,4*(BeamSegment+1)+2+4,:))+node.coordinates_mat(4*(BeamSegment+1)+2+4,:)';
    vector(i,:)=node1-node2; 
    angleTwist(i)=atan(vector(i,3)/vector(i,2))*180/pi;
end

frequency=(f2-f1)*time(1:CorpFactor:end)/totalSweepTime+f1;

figure
plot(frequency,angleTwist)
xlabel('frequency (Hz)');
ylabel('Twisting angle');
title('Fast Axis Bent Mirror')


%% Check bending also
angleBend=zeros(step/CorpFactor,1);
vector=zeros(step/CorpFactor,3);

for i=1:step/CorpFactor
    node1=squeeze(UhisCorp(i,4*(BeamSegment+1)+2+7,:))+node.coordinates_mat(4*(BeamSegment+1)+2+7,:)';
    node2=squeeze(UhisCorp(i,4*(BeamSegment+1)+2+2,:))+node.coordinates_mat(4*(BeamSegment+1)+2+2,:)';
    vector(i,:)=node1-node2;    
    angleBend(i)=sign(vector(i,1))*atan(vector(i,3)/vector(i,1))*180/pi;
    if vector(i,1)<0
        angleBend(i)=180-angleBend(i);
    end

end

figure
plot(frequency,angleBend)
xlabel('Frequency (Hz)');
ylabel('Bending Angle');
title('Bending Motoin of Micro Mirror')

UhisPlot=UhisCorp(1:50:end,:,:);
plots.fileName='MicroMirror.gif';
plots.displayRange=1.5*10^-3;
plots.displayRangeRatio=0.6;
plots.viewAngle1=10;
plots.viewAngle2=20;
plots.Plot_DeformedHis(UhisPlot);






