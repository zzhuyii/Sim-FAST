clear all;
clc;
close all;

%% Define Geometry, HSS 4X3X5/16 A500 Grade C Fy=50ksi
H=2;
HA=2;
W=2;
L=2;
l=0.3;

N=8;

barA=0.0023; % 3.52 in^2
barE=2*10^11;
panel_E=2*10^11;
panel_t=0.05;
panel_v=0.3;
activeBarE=2*10^11; % 80*10^9;

% I=1/12*0.01^4;
% kspr=3*barE*I/L2*100000;

node=Elements_Nodes;
bar=Vec_Elements_Bars;
actBar=CD_Elements_Bars;
cst=Vec_Elements_CST;

assembly=Assembly_Truss_Rolling_Bridge;
assembly.node=node;
assembly.bar=bar;
assembly.actBar=actBar;
assembly.cst=cst;

%% Define the nodal coordinates
node.coordinates_mat=[];
node.coordinates_mat=[node.coordinates_mat;
    0    0   0;
    L/2  0   H;
    L    0   0;
    0    W   0;
    L    W   0;
    L/2  W   H;
    L    0   H;
    L    W   H;
    ];

for i=2:N-1
    node.coordinates_mat=[node.coordinates_mat;
        L*i      0   0;
        L*i-L/2  0   H;
        L*i      W   0;
        L*i-L/2  W   H;
        L*i      0   H;
        L*i      W   H;
        ];
end

node.coordinates_mat=[node.coordinates_mat;
    L*N      0   0;
    L*N-L/2  0   H;
    L*N      W   0;
    L*N-L/2  W   H;
    ];

% Set up the plotting function for inspection
plots=Plot_Truss_Rolling_Bridge();
plots.assembly=assembly;

% We will plot for the Rolling Bridge
plots.displayRange=[-0.5;2*N+0.5;-0.5;2.5;-0.5;2.5]; 

plots.viewAngle1=20;
plots.viewAngle2=20;


% Plot the nodal coordinates for inspection
plots.Plot_Shape_Node_Number()


%% Define how panels are designed
cst.node_ijk_mat=[cst.node_ijk_mat;
    1  3  4;
    3  4  5;
    ];

for i=2:N
    cst.node_ijk_mat=[cst.node_ijk_mat;
        3+(i-2)*6  5+(i-2)*6  9+(i-2)*6;
        5+(i-2)*6  9+(i-2)*6  11+(i-2)*6;
        ];
end    

cstNum=size(cst.node_ijk_mat,1);
cst.v_vec=panel_v*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.t_vec=panel_t*ones(cstNum,1);

plots.Plot_Shape_CST_Number();


%% Define how normal bars are connected
% First we design the normal bar
bar.node_ij_mat=[bar.node_ij_mat;
    1, 2;
    1, 3;
    2, 3;
    4, 5;
    4, 6;
    5, 6;
    2, 7;
    6, 8;
    1, 4;
    3, 5;
    3, 4;
    ];

for i=2:N-1
    bar.node_ij_mat=[bar.node_ij_mat;
        3+(i-2)*6,  9+(i-2)*6;
        3+(i-2)*6,  10+(i-2)*6;
        9+(i-2)*6,  10+(i-2)*6;
        5+(i-2)*6,  11+(i-2)*6;
        5+(i-2)*6,  12+(i-2)*6;
        11+(i-2)*6, 12+(i-2)*6;
        7+(i-2)*6,  10+(i-2)*6;
        8+(i-2)*6,  12+(i-2)*6;
        10+(i-2)*6, 13+(i-2)*6;
        12+(i-2)*6, 14+(i-2)*6;
        9+(i-2)*6,  11+(i-2)*6;
        5+(i-2)*6,  9+(i-2)*6;
        ];
end

bar.node_ij_mat=[bar.node_ij_mat;
    3+(N-2)*6,  9+(N-2)*6;
    3+(N-2)*6,  10+(N-2)*6;
    9+(N-2)*6,  10+(N-2)*6;
    5+(N-2)*6,  11+(N-2)*6;
    5+(N-2)*6,  12+(N-2)*6;
    11+(N-2)*6, 12+(N-2)*6;
    7+(N-2)*6,  10+(N-2)*6;
    8+(N-2)*6,  12+(N-2)*6;
    5+(N-2)*6,  9+(N-2)*6;
    9+(N-2)*6,  11+(N-2)*6;
    ];

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);

plots.Plot_Shape_Bar_Number();

%% Define how actuator bars are connected
% Next we design the active bars
for i=1:N-1
    actBar.node_ij_mat=[actBar.node_ij_mat;
        3+(i-1)*6, 7+(i-1)*6;
        5+(i-1)*6, 8+(i-1)*6;
        ];
end

actBarNum=size(actBar.node_ij_mat,1);
actBar.A_vec=barA*ones(actBarNum,1);
actBar.E_vec=activeBarE*ones(actBarNum,1);

plots.Plot_Shape_ActBar_Number();

% Initialize the entire assembly 
assembly.Initialize_Assembly();

%% Calculate Self-weight of the Bridge
% Steel properties
rho_steel=7850;         % Density of steel (kg/m^3)
g=9.81;                 % Gravitational acceleration (m/s^2)

% Bars (including normal bars and actuator bars)
% Bar cross-section area
A_bar=barA;             % m^2, should match bar.A_vec

% 1. Total length of all normal bar elements
L_bar_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;
for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    p1=coords(n1,:);
    p2=coords(n2,:);
    len=norm(p1-p2);
    L_bar_total=L_bar_total+len;
end

% 2. Total length of all actuator bar elements
L_actbar_total=0;
if isfield(actBar,'node_ij_mat')
    actBarNodeMat=actBar.node_ij_mat;
    for i=1:size(actBarNodeMat,1)
        n1=actBarNodeMat(i,1);
        n2=actBarNodeMat(i,2);
        p1=coords(n1,:);
        p2=coords(n2,:);
        len=norm(p1-p2);
        L_actbar_total=L_actbar_total+len;
    end
end

% 3. Total bar length (all bars)
L_total=L_bar_total+L_actbar_total;

% 4. Total weight of all bars (Newtons)
W_bar=A_bar*L_total*rho_steel*g;

% CST panels
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
    area=sqrt(s*(s-a)*(s-b)*(s-c));
    A_cst_total=A_cst_total+area;
end

t_cst=panel_t;  % Thickness of CST panel (m), matches cst.t_vec
% Weight of CST panels (Newtons)
W_cst=A_cst_total*t_cst*rho_steel*g;

% ----- Total self-weight -----
W_total=W_bar+W_cst;

% Output results
fprintf('-----------------------------\n');
fprintf('Total length of all bars: %.2f m\n', L_total);
fprintf('Total area of all CST panels: %.2f m^2\n', A_cst_total);
fprintf('Total bar weight: %.2f N\n', W_bar);
fprintf('Total CST panel weight: %.2f N\n', W_cst);
fprintf('Total self-weight of the bridge: %.2f N\n', W_total);
fprintf('-----------------------------\n');



%% Distributed load along full length on bridge bottom
nr=Solver_NR_Loading;
nr.assembly=assembly;

%% Supports 
coords=node.coordinates_mat;
tol=1e-9;

nodeNum=size(coords,1);
nodeNumVec=(1:nodeNum)';

% identify four bottom corners by coordinates
xmin=min(coords(:,1));
xmax=max(coords(:,1));

isBottom=abs(coords(:,3)-0)<tol;
isY0=abs(coords(:,2)-0)<tol;
isYW=abs(coords(:,2)-W)<tol;
isXmin=abs(coords(:,1)-xmin)<tol;
isXmax=abs(coords(:,1)-xmax)<tol;

A=find(isBottom & isXmin & isY0, 1);   % left-front-bottom corner
B=find(isBottom & isXmin & isYW, 1);   % left-back-bottom corner
C=find(isBottom & isXmax & isY0, 1);   % right-front-bottom corner
D=find(isBottom & isXmax & isYW, 1);   % right-back-bottom corner

if any(cellfun(@isempty, {A,B,C,D}))
    error('Failed to find four bottom corner nodes. Check geometry/tol/W.');
end

fprintf('Corner supports: A=%d, B=%d, C=%d, D=%d\n', A, B, C, D);

% Stabilization like your original:
% lock uy for all nodes (3 dof per node: [ux uy uz] => [0 1 0])
nr.supp = [nodeNumVec, zeros(nodeNum,1), ones(nodeNum,1), zeros(nodeNum,1)];

% fully fix the 4 corners (ux=uy=uz=1)
nr.supp(A,2:4)=1;
nr.supp(B,2:4)=1;
nr.supp(C,2:4)=1;
nr.supp(D,2:4)=1;


force=1000;   % N  (IMPORTANT: interpreted as "total load added PER increment")
step=50;     % number of increments (final total load = force*step)

%  1) Bridge length along x 
L_total=L*N;              % m
q=force/L_total;       % N/m  (line load PER increment, consistent with "force per increment")

% 2) Find bottom-edge nodes (z=0 and y=0 or y=W)
% coords=node.coordinates_mat;
% tol=1e-9;

isBottomZ=abs(coords(:,3)-0)<tol;
isEdgeY0=abs(coords(:,2)-0)<tol;
isEdgeYW=abs(coords(:,2)-W)<tol;

bottom=find(isBottomZ&(isEdgeY0|isEdgeYW));   % bottom edge nodes (full length)

if numel(bottom)<2
    error('Bottom node set is too small. Check y/z criteria or geometry.');
end

% 3) Group bottom nodes by x-section (robust to two-edge nodes at same x)
x_all=coords(bottom,1);

tolX=1e-9;
x_round=round(x_all/tolX)*tolX;

[xu, ~, g]=unique(x_round, 'stable');   % xu: unique x positions (not necessarily sorted)
[xu, ord]=sort(xu);                    % sort x sections
g=ord(g);                      % remap group IDs after sorting xu

nsec=numel(xu);
if nsec<2
    error('Not enough x-sections to compute tributary lengths.');
end

% 4) Tributary length per x-section
ell_sec=zeros(nsec,1);
for s=1:nsec
    if s==1
        ell_sec(s)=(xu(s+1)-xu(s))/2;
    elseif s==nsec
        ell_sec(s)=(xu(s)-xu(s-1))/2;
    else
        ell_sec(s)=(xu(s+1)-xu(s-1))/2;
    end
end

% --- 5) Equivalent nodal forces PER increment (sum over all bottom nodes = -force) ---
Fz_base=zeros(numel(bottom),1);
for s=1:nsec
    ids=find(g==s);          % indices within 'bottom' at same x-section
    Fsec=-q*ell_sec(s);       % total resultant at this x-section (PER increment)
    Fz_base(ids)=Fsec/numel(ids);   % split among nodes at same x (usually 2: y=0 and y=W)
end

fprintf('Check sum(Fz_base) = %.6f N (target = %.6f N)\n', sum(Fz_base), -force);

% 6) Assign load table: [nodeID, Fx, Fy, Fz] 
nr.load=[bottom(:), zeros(numel(bottom),2), Fz_base(:)];

% 7) Nonlinear solver settings
nr.increStep=step;
nr.iterMax =50;
nr.tol=1e-5;


% 8) Solve
Uhis=nr.Solve;



%  0) Final total external load 
P_final=force*step;    % N (final total vertical load)

% Optional sanity checks
fprintf('Postproc: Per-increment load = %.6f N, steps = %d, P_final = %.6f N\n', ...
        force, step, P_final);
fprintf('Postproc: Sum of nr.load(:,4) = %.6f N (should be about -force per increment)\n', ...
        sum(nr.load(:,4)));

% 1) Final displacement state
U_end=squeeze(Uhis(end,:,:));   % [nodeNum x 3]

% 2) Bar strain and internal axial force
truss_strain=bar.Solve_Strain(node, U_end);                         % strain in each bar
internal_force=truss_strain.*(bar.E_vec).*(bar.A_vec);            % axial force in each bar (N)

% 3) Max internal force
[maxBarForce, idxMax]=max(abs(internal_force));
fprintf('Max |bar force|=%.6f N at bar #%d\n', maxBarForce, idxMax);

% 4) Ultimate capacity (simple tensile stress limit)
sigma_u=300e6;          % Pa
barFailureForce=sigma_u*barA; % N

% 5) Total bar length 
barLtotal=sum(bar.L0_vec);

% 6) Stiffness evaluation (midspan bottom nodes)
x_mid =0.5*L_total;

% nearest x-section
[~, sMid]=min(abs(xu-x_mid));
mid_ids_local = find(g==sMid);        % indices within 'bottom' at that x-section
mid_nodes=bottom(mid_ids_local);      % actual node IDs at midspan bottom

% average downward displacement at those nodes
Uaverage=-mean(U_end(mid_nodes,3));   % positive downward

% stiffness based on final total load
Kstiff=P_final/Uaverage;            % N/m

fprintf('Midspan nodes used for stiffness: %s\n', mat2str(mid_nodes(:)'));
fprintf('Uaverage = %.6e m, Kstiff = %.6e N/m\n', Uaverage, Kstiff);

% 7) Stress-like visualization field
% This keeps your original plotting logic: scale so the critical bar reaches barFailureForce
bar_stress=truss_strain .* (bar.E_vec) * (barFailureForce / maxBarForce);
plots.Plot_Shape_Bar_Stress(bar_stress);

% 8) Failure load prediction by linear scaling 
% Assuming internal forces scale linearly with load:
loadatfail=P_final * (barFailureForce / maxBarForce);   % N

fprintf('Failure load is %.3f kN\n', loadatfail/1000);
fprintf('Total bar length is %.3f m\n', barLtotal);
fprintf('Stiffness is %.3e N/m\n', Kstiff);

% 9) Optional efficiency (if W_bar exists)
if exist('W_bar','var') && ~isempty(W_bar) && W_bar>0
    loadEff = loadatfail / W_bar;
    fprintf('Load efficiency (loadatfail/W_bar) = %.6f\n', loadEff);
end

% 10) Optional: report top critical bars
nb=numel(internal_force);
topk=min(10, nb);
[~, ord]=sort(abs(internal_force), 'descend');

fprintf('Worst %d bars by |axial force|:\n', topk);
for ii=1:topk
    b=ord(ii);
    fprintf('#%d: Bar %d, N = %+ .3f kN\n', ii, b, internal_force(b)/1e3);
end

% 11) Optional: compute actual axial stress in each bar (Pa)
bar_sigma=internal_force ./ bar.A_vec;   % Pa
sigma_max=max(abs(bar_sigma));
fprintf('Max |axial stress| = %.3f MPa\n', sigma_max/1e6);




%% Set up the self actuation solver
nr=Solver_NR_Loading;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

[T,K]=assembly.Solve_FK(zeros(size(node.coordinates_mat)));
spy(K)

nr.assembly=assembly;
nr.supp=[nodeNumVec,zeros(nodeNum,1),ones(nodeNum,1),zeros(nodeNum,1)];
nr.supp(1,2:4)=ones(1,3);
nr.supp(4,2:4)=ones(1,3);
nr.supp(6*N-3,2:4)=ones(1,3); % 21  45
nr.supp(6*N-1,2:4)=ones(1,3); % 23  47

P_total=10000;
step=10;

nr.load=[];
for i=1:N-1
    nr.load=[nr.load;
        3+(i-1)*6 0 0 -P_total/2/(N-1);
        5+(i-1)*6 0 0 -P_total/2/(N-1);
        ];
end

nr.increStep=step;
nr.iterMax=30;
nr.tol=10^-5;

Uhis=nr.Solve();

% Plot the deformed shape
plots.Plot_Deformed_Shape(20*squeeze(Uhis(end,:,:)));

% Also plot the deformation history
plots.fileName="Rolling_Bridge_Load.gif";
plots.Plot_Deformed_His(Uhis);

% Failure load computation
truss_strain=bar.Solve_Strain(node,squeeze(Uhis(end,:,:))*50);
internal_force=(truss_strain).*(bar.E_vec).*(bar.A_vec);

% Find the maximum bar internal force
[maxBarForce,index]=max(abs(internal_force));

% Find failure force for the bar
sigma_u=300*10^6;
barFailureForce=sigma_u*barA;

% Find total bar length
barLtotal=sum(bar.L0_vec);

% Find Stiffness
Uaverage=-mean(squeeze(Uhis(end,[3*N-3,3*N-1],3)));
Kstiff=step*force/Uaverage;

% Plot failure stress
bar_stress=(truss_strain).*(bar.E_vec)*barFailureForce/maxBarForce;
plots.Plot_Shape_Bar_Stress(bar_stress)

% Find the relationship betweeen the bar internal forces and load
loadatfail=force*step*barFailureForce/maxBarForce;
fprintf('Failure load is %d kN \n',  loadatfail/1000);
fprintf('Total bar length is %d m \n',  barLtotal);
fprintf('Stiffness is %d N/m \n',  Kstiff);
fprintf('Normal bars: %d\n', barNum);
fprintf('Actuator bars: %d\n', actBarNum);
fprintf('Total bars: %d\n', barNum + actBarNum);



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
K =1.0;                
Lc=K.*L0_vec;

% 2) r, r = 0.5*sqrt(A/pi)
r=0.5*sqrt(A./pi);
r=max(r,1e-9); % prevent division by zero

% 3) yield stress
Fy=250e6;             % Pa（need to be changed）

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