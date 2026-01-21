clear all
close all
clc
tic

%% Define Geometry, HSS 4X3X5/16 A500 Grade C Fy=50ksi
N=8; 
L=2; 
W=2;
H=2;
gap=0;

barA=0.0023; % 3.52 in^2
barE=2*10^11;
panel_E=2*10^11;
panel_t=0.05;
panel_v=0.3;


node=Elements_Nodes;
node.coordinates_mat=[node.coordinates_mat;
    0, 0, 0; % 1
    0, W, 0; % 2
    0, 0, H; % 3
    0, W, H; % 4
    ]; 

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        (L)*(i-1)+L/2, 0, 0;
        (L)*(i-1)+L/2, 0, gap;
        (L)*(i-1)+L/2, W, 0;
        (L)*(i-1)+L/2, W, gap;

        (L)*(i-1)+L/2, 0, H;
        (L)*(i-1)+L/2, gap, H;
        (L)*(i-1)+L/2, W, H;
        (L)*(i-1)+L/2, W, H-gap;

        (L)*(i-1)+L/2, W/2, 0;
        (L)*(i-1)+L/2, W/2, H;
        (L)*(i-1)+L/2, 0, H/2;
        (L)*(i-1)+L/2, W, H/2;

        (L)*(i-1)+L, 0, 0;
        (L)*(i-1)+L, W, 0;
        (L)*(i-1)+L, 0, H;
        (L)*(i-1)+L, W, H;
        ];
end


%% Define assembly
assembly=Assembly_Kirigami_Truss;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;
bar=Vec_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_4N=rot_spr_4N;

%% Define Plotting Functions
plots=Plot_Kirigami_Truss;
plots.assembly=assembly;
plots.displayRange=[-0.5; 2*N+0.5; -0.5; 2.5; -0.5; 2.5];

plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_Node_Number;


%% Define Triangle
for i=1:N
     cst.node_ijk_mat=[cst.node_ijk_mat;
        16*(i-1)+1    16*(i-1)+2    16*(i-1)+13;
        16*(i-1)+1    16*(i-1)+5    16*(i-1)+13;
        16*(i-1)+2    16*(i-1)+7    16*(i-1)+13;
        16*(i-1)+13    16*(i-1)+17    16*(i-1)+18;
        16*(i-1)+6    16*(i-1)+13    16*(i-1)+17;
        16*(i-1)+8    16*(i-1)+13    16*(i-1)+18;
        ];
end

cstNum=size(cst.node_ijk_mat,1);
cst.t_vec=panel_t*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.v_vec=panel_v*ones(cstNum,1);

plots.Plot_Shape_CST_Number;

%% Define bar
bar.node_ij_mat=[bar.node_ij_mat;
    1 2;
    1 3;
    3 4;
    2 4;
    ];

for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        16*(i-1)+3   16*(i-1)+15;
        16*(i-1)+1  16*(i-1)+15;
        16*(i-1)+1   16*(i-1)+5;
        16*(i-1)+5   16*(i-1)+15;
        16*(i-1)+3   16*(i-1)+9;
        16*(i-1)+9  16*(i-1)+15;
        16*(i-1)+10   16*(i-1)+15;
        16*(i-1)+10   16*(i-1)+19;
        16*(i-1)+15   16*(i-1)+19;
        16*(i-1)+15  16*(i-1)+17;
        16*(i-1)+6   16*(i-1)+15;
        16*(i-1)+6   16*(i-1)+17;

        16*(i-1)+3   16*(i-1)+9;
        16*(i-1)+3  16*(i-1)+14;
        16*(i-1)+9   16*(i-1)+14;
        16*(i-1)+4   16*(i-1)+14;
        16*(i-1)+4   16*(i-1)+11;
        16*(i-1)+11  16*(i-1)+14;
        16*(i-1)+10   16*(i-1)+14;
        16*(i-1)+10   16*(i-1)+19;
        16*(i-1)+14   16*(i-1)+19;
        16*(i-1)+14  16*(i-1)+20;
        16*(i-1)+12   16*(i-1)+14;
        16*(i-1)+12   16*(i-1)+20;

        16*(i-1)+2   16*(i-1)+16;
        16*(i-1)+2  16*(i-1)+7;
        16*(i-1)+7   16*(i-1)+16;
        16*(i-1)+4   16*(i-1)+16;
        16*(i-1)+4   16*(i-1)+11;
        16*(i-1)+11  16*(i-1)+16;
        16*(i-1)+12   16*(i-1)+16;
        16*(i-1)+12   16*(i-1)+20;
        16*(i-1)+16   16*(i-1)+20;
        16*(i-1)+16  16*(i-1)+18;
        16*(i-1)+8   16*(i-1)+16;
        16*(i-1)+8   16*(i-1)+18;

        16*(i-1)+1   16*(i-1)+5;
        16*(i-1)+1  16*(i-1)+13;
        16*(i-1)+5   16*(i-1)+13;
        16*(i-1)+2   16*(i-1)+13;
        16*(i-1)+2   16*(i-1)+7;
        16*(i-1)+7  16*(i-1)+13;
        16*(i-1)+8   16*(i-1)+13;
        16*(i-1)+8   16*(i-1)+18;
        16*(i-1)+13   16*(i-1)+18;
        16*(i-1)+13  16*(i-1)+17;
        16*(i-1)+6   16*(i-1)+13;
        16*(i-1)+6   16*(i-1)+17;

        16*(i-1)+17   16*(i-1)+19;
        16*(i-1)+17  16*(i-1)+18;
        16*(i-1)+19   16*(i-1)+20;
        16*(i-1)+18   16*(i-1)+20;
        ];
end

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);
plots.Plot_Shape_Bar_Number();

%% Define Rotational Spring
for i=1:N    
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        16*(i-1)+1  16*(i-1)+3  16*(i-1)+15  16*(i-1)+9;
        16*(i-1)+3  16*(i-1)+1  16*(i-1)+15  16*(i-1)+5;
        16*(i-1)+17  16*(i-1)+19  16*(i-1)+15  16*(i-1)+10;
        16*(i-1)+19  16*(i-1)+17  16*(i-1)+15  16*(i-1)+6;
        16*(i-1)+3  16*(i-1)+9  16*(i-1)+15  16*(i-1)+19;
        16*(i-1)+1  16*(i-1)+5  16*(i-1)+15  16*(i-1)+17;
        16*(i-1)+17  16*(i-1)+6  16*(i-1)+18  16*(i-1)+1;
        16*(i-1)+19  16*(i-1)+10  16*(i-1)+15  16*(i-1)+3; %8

        16*(i-1)+4  16*(i-1)+3  16*(i-1)+14  16*(i-1)+9;
        16*(i-1)+3  16*(i-1)+14  16*(i-1)+9  16*(i-1)+19;
        16*(i-1)+3  16*(i-1)+4  16*(i-1)+14  16*(i-1)+11;
        16*(i-1)+4  16*(i-1)+11  16*(i-1)+14  16*(i-1)+20;
        16*(i-1)+20  16*(i-1)+19  16*(i-1)+14  16*(i-1)+10;
        16*(i-1)+19  16*(i-1)+14  16*(i-1)+10  16*(i-1)+3;
        16*(i-1)+19  16*(i-1)+20  16*(i-1)+14  16*(i-1)+12;
        16*(i-1)+20  16*(i-1)+12  16*(i-1)+14  16*(i-1)+4;

        16*(i-1)+2  16*(i-1)+4  16*(i-1)+16  16*(i-1)+11;
        16*(i-1)+4  16*(i-1)+11  16*(i-1)+16  16*(i-1)+20;
        16*(i-1)+4  16*(i-1)+2  16*(i-1)+16  16*(i-1)+7;
        16*(i-1)+2  16*(i-1)+16  16*(i-1)+7  16*(i-1)+18;
        16*(i-1)+18  16*(i-1)+20  16*(i-1)+16  16*(i-1)+12;
        16*(i-1)+20  16*(i-1)+12  16*(i-1)+16  16*(i-1)+4;
        16*(i-1)+20  16*(i-1)+18  16*(i-1)+16  16*(i-1)+8;
        16*(i-1)+18  16*(i-1)+16  16*(i-1)+8  16*(i-1)+2;

        16*(i-1)+2  16*(i-1)+1  16*(i-1)+13  16*(i-1)+5;
        16*(i-1)+1  16*(i-1)+13  16*(i-1)+5  16*(i-1)+17;
        16*(i-1)+1  16*(i-1)+2  16*(i-1)+13  16*(i-1)+7;
        16*(i-1)+2  16*(i-1)+7  16*(i-1)+13  16*(i-1)+18;
        16*(i-1)+17  16*(i-1)+18  16*(i-1)+13  16*(i-1)+8;
        16*(i-1)+18  16*(i-1)+8  16*(i-1)+13  16*(i-1)+2;
        16*(i-1)+18  16*(i-1)+17  16*(i-1)+13  16*(i-1)+6;
        16*(i-1)+17  16*(i-1)+6  16*(i-1)+13  16*(i-1)+1;
        ];
        
end


rotNum=size(rot_spr_4N.node_ijkl_mat);
rotNum=rotNum(1);

rot_spr_4N.rot_spr_K_vec=ones(rotNum,1);

% factor=100;
% for i=1:N+1
%     rot_spr_4N.rot_spr_K_vec((i-1)*40+2)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+2);
%     rot_spr_4N.rot_spr_K_vec((i-1)*40+4)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+4);
%     rot_spr_4N.rot_spr_K_vec((i-1)*40+6)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+6);
%     rot_spr_4N.rot_spr_K_vec((i-1)*40+8)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+8);
% end

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;

assembly.Initialize_Assembly;


%% Calculate self-weight 
rho_steel=7850;         % kg/m^3
A_bar=barA;             % m^2, (same as bar.A_vec)

% a. Find total length of all bar elements
L_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;

for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    len=norm(coords(n1,:)-coords(n2,:));
    L_total=L_total+len;
end

W_bar=A_bar*L_total*rho_steel*9.81;   % unit:N

% c. Find total area of all CST elements
A_cst_total=0;
cstNodeMat=cst.node_ijk_mat;
for i=1:size(cstNodeMat,1)
    n1=cstNodeMat(i,1);
    n2=cstNodeMat(i,2);
    n3=cstNodeMat(i,3);
    p1=coords(n1,:);
    p2=coords(n2,:);
    p3=coords(n3,:);
    % Heron's formula for triangle area
    a=norm(p1-p2);
    b=norm(p2-p3);
    c=norm(p3-p1);
    s=(a+b+c)/2;
    area = sqrt(s*(s-a)*(s-b)*(s-c));
    A_cst_total = A_cst_total + area;
end

t_cst=panel_t; 
W_cst=A_cst_total*t_cst*rho_steel*9.81;   % unit:N

% self-weight
W_total=W_bar+W_cst;

fprintf('Total bar length: %.2f m\n', L_total);
fprintf('Total CST area: %.2f m^2\n', A_cst_total);
fprintf('Bar weight: %.2f N\n', W_bar);
fprintf('CST weight: %.2f N\n', W_cst);
fprintf('Total self-weight: %.2f N\n', W_total);


%% Set up solver + Distributed load on bottom nodes (full-length)
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[1 1 1 1;
         2 1 1 1;
         16*N+1 1 1 1;
         16*N+2 1 1 1;];

P_total=10000;  % N (total vertical load)
step=10;

nr.load=[];
for i=1:N-1
    nr.load=[nr.load;
        5+(i-1)*16 0 0 -P_total/2/(N-1);
        7+(i-1)*16 0 0 -P_total/2/(N-1);
        17+(i-1)*16 0 0 -P_total/2/(N-1);
        18+(i-1)*16 0 0 -P_total/2/(N-1);
        ];
end

nr.increStep=1;
nr.iterMax=20;
nr.tol=1*10^-5;

Uhis=nr.Solve;
          

plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:))*50)


% Post-processing: bar strain, internal force, failure load, stiffness
truss_strain=bar.Solve_Strain(node,squeeze(Uhis(end,:,:)));
internal_force=(truss_strain).*(bar.E_vec).*(bar.A_vec);


% Find the maximum bar internal force
[maxBarForce,index]=max(abs(internal_force));


% Find failure force for the bar
sigma_u=300*10^6; % ultimate stress, 300 MPa
barFailureForce=sigma_u*barA; % N


% Find total bar length
barLtotal=sum(bar.L0_vec);

% Find Stiffness
Uaverage=-mean(squeeze(Uhis(end,[N*8+1,N*8+2],3)));
Kstiff=step*P_total/Uaverage;


% Plot failure stress
bar_stress=(truss_strain).*(bar.E_vec);
plots.Plot_Shape_Bar_Stress(bar_stress)


% Find the relationship betweeen the bar internal forces and load
loadatfail=P_total*step*barFailureForce/maxBarForce;
fprintf('Failure load is %d kN \n',  loadatfail/1000);
fprintf('Total bar length is %d m \n',  barLtotal);
fprintf('Stiffness is %d N/m \n',  Kstiff);
fprintf('Total bars: %d\n', barNum);

loadEff=loadatfail/W_bar;
fprintf('Load efficiency (FailureLoad/SelfWeight) = %.3f \n', loadEff);



%% Evaluate Member
AxialForce=internal_force(:); % N
A=bar.A_vec(:); % m^2
E=bar.E_vec(:); % Pa
nb=numel(A);

% 1) effective length KL
if ~isfield(bar,'L0_vec')||isempty(bar.L0_vec)
    L0_vec=zeros(nb,1);
    for k=1:nb
        n1=bar.node_ij_mat(k,1);
        n2=bar.node_ij_mat(k,2);
        L0_vec(k)=norm(assembly.node.coordinates_mat(n1,:) - ...
                         assembly.node.coordinates_mat(n2,:));
    end
else
    L0_vec=bar.L0_vec(:);
end
K=1.0;                
Lc=K.*L0_vec;

% 2) r, r = 0.5*sqrt(A/pi)
r=0.5*sqrt(A./pi);
r=max(r,1e-9); % prevent division by zero

% 3) yield stress
Fy=345e6;             % Pa（50ksi）

% 4) evaluate
passYN=false(nb,1);
util=NaN(nb,1);
modeStr=cell(nb,1);
Pn=NaN(nb,1);
slender=NaN(nb,1);   % KL/r
Fe=NaN(nb,1);   % Euler stress (Pa)
Fcr=NaN(nb,1);   % Critical stress per AISC (Pa)

tiny=1e-12;
for i=1:nb
    Ni=AxialForce(i);
    Ai=A(i);
    Ei=E(i);
    Lci=Lc(i);
    ri=r(i);

    if Ni>0
        % in tension 
        Pn_i=Fy*Ai;
        modeStr{i}='Tension-Yield';
        slender(i)=NaN;  
        Fe(i)=NaN; 
        Fcr(i)=NaN;   % don't calculate buckling in tension
    else
        % in compression
        slender=Lci/ri;
        Fe=(pi^2*Ei)/(slender^2);
        lambda_lim=4.71*sqrt(Ei/Fy);
        if slender<=lambda_lim
            Fcr=(0.658)^(Fy/Fe)*Fy;   % inelastic buckling
        else
            Fcr=0.877*Fe;             % elastic bukling
        end
        Pn_i=Fcr*Ai;
        modeStr{i}='Compression-Buckling';
    end

    Pn(i)=Pn_i;
    util(i)=abs(Ni)/max(Pn_i,tiny);
    passYN(i)=util(i)<=1.0;
end

% 5) print
fprintf('Bars passed: %d / %d\n', sum(passYN), nb);

[util_sorted, idx]=sort(util, 'descend');
topk=min(10, nb);
fprintf('Worst %d bars (by utilization):\n', topk);
for ii=1:topk
    b=idx(ii);
    fprintf('#%d: N=%.2f kN, util=%.3f, mode=%s, Pn=%.2f kN\n', ...
        b, AxialForce(b)/1e3, util_sorted(ii), modeStr{b}, Pn(b)/1e3);
end