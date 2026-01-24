clear all
close all
clc
tic

%% Define Geometry
N=4;
L=0.875; %0.875(N=4) 0.8875(N=8)
W=1;
H=1;
w=0.1;
gap=0;


barA=0.00227; % 0.0004; % m^2
barE=200*10^9; % Pa
panel_E=200*10^9;
panel_t=0.05;
panel_v=0.2;

node=Elements_Nodes;
node.coordinates_mat=[node.coordinates_mat;
        -w, 0, 0; % 1
        -w, W, 0; % 2
        -w, 0, H; % 3
        -w, W, H;]; % 4

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        (w+L)*(i-1), 0, 0;
        (w+L)*(i-1), W, 0;
        (w+L)*(i-1), 0, H;
        (w+L)*(i-1), W, H;

        (w+L)*(i-1)+L/2, 0, 0;
        (w+L)*(i-1)+L/2, 0, gap;
        (w+L)*(i-1)+L/2, W, 0;
        (w+L)*(i-1)+L/2, W, gap;

        (w+L)*(i-1)+L/2, 0, H;
        (w+L)*(i-1)+L/2, gap, H;
        (w+L)*(i-1)+L/2, W, H;
        (w+L)*(i-1)+L/2, W, H-gap;

        (w+L)*(i-1)+L/2, W/2, 0;
        (w+L)*(i-1)+L/2, W/2, H;
        (w+L)*(i-1)+L/2, 0, H/2;
        (w+L)*(i-1)+L/2, W, H/2;

        (w+L)*(i-1)+L, 0, 0;
        (w+L)*(i-1)+L, W, 0;
        (w+L)*(i-1)+L, 0, H;
        (w+L)*(i-1)+L, W, H;];
end

node.coordinates_mat=[node.coordinates_mat;
        (w+L)*N, 0, 0;
        (w+L)*N, W, 0;
        (w+L)*N, 0, H;
        (w+L)*N, W, H;
        ];


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
plots.displayRange=[-0.5; N+0.5; -0.5; 1.5; -0.5; 1.5];

plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_Node_Number;


%% Define Triangle
for i=1:N+1
     cst.node_ijk_mat=[cst.node_ijk_mat;
        20*(i-1)+1    20*(i-1)+3    20*(i-1)+7;
        20*(i-1)+1    20*(i-1)+5    20*(i-1)+7;
        20*(i-1)+3    20*(i-1)+4    20*(i-1)+8;
        20*(i-1)+3    20*(i-1)+7    20*(i-1)+8;
        20*(i-1)+2    20*(i-1)+4    20*(i-1)+6;
        20*(i-1)+4    20*(i-1)+6    20*(i-1)+8;
        20*(i-1)+3    20*(i-1)+4    20*(i-1)+8;
        20*(i-1)+1    20*(i-1)+2    20*(i-1)+6;
        20*(i-1)+1    20*(i-1)+5    20*(i-1)+6;];
end

for i=1:N
     cst.node_ijk_mat=[cst.node_ijk_mat;
        20*(i-1)+5    20*(i-1)+6    20*(i-1)+17;
        20*(i-1)+5    20*(i-1)+9    20*(i-1)+17;
        20*(i-1)+11    20*(i-1)+6    20*(i-1)+17;
        20*(i-1)+9    20*(i-1)+21    20*(i-1)+17;
        20*(i-1)+11    20*(i-1)+22    20*(i-1)+17;
        20*(i-1)+21    20*(i-1)+22    20*(i-1)+17;];
end

cstNum=size(cst.node_ijk_mat,1);
cst.t_vec=panel_t*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.v_vec=panel_v*ones(cstNum,1);

plots.Plot_Shape_CST_Number;

%% Define bar
for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        20*(i-1)+5   20*(i-1)+10;
        20*(i-1)+19  20*(i-1)+10;
        20*(i-1)+5   20*(i-1)+19;
        20*(i-1)+7   20*(i-1)+13;
        20*(i-1)+7   20*(i-1)+19;
        20*(i-1)+13   20*(i-1)+19;
        20*(i-1)+13   20*(i-1)+23;
        20*(i-1)+19   20*(i-1)+23;
        20*(i-1)+19   20*(i-1)+21;
        20*(i-1)+10   20*(i-1)+21; %10

        20*(i-1)+7   20*(i-1)+14;
        20*(i-1)+7   20*(i-1)+18;
        20*(i-1)+14   20*(i-1)+18;
        20*(i-1)+8   20*(i-1)+15;
        20*(i-1)+15   20*(i-1)+18;
        20*(i-1)+8   20*(i-1)+18;
        20*(i-1)+14   20*(i-1)+23;
        20*(i-1)+18   20*(i-1)+23;
        20*(i-1)+18   20*(i-1)+24;
        20*(i-1)+15   20*(i-1)+24; %20

        20*(i-1)+8   20*(i-1)+16;
        20*(i-1)+16   20*(i-1)+20;
        20*(i-1)+8   20*(i-1)+20;
        20*(i-1)+24   20*(i-1)+16;
        20*(i-1)+24   20*(i-1)+20;
        20*(i-1)+6   20*(i-1)+20;
        20*(i-1)+6   20*(i-1)+12;
        20*(i-1)+12   20*(i-1)+20;
        20*(i-1)+22   20*(i-1)+12;
        20*(i-1)+22   20*(i-1)+20; %30

        20*(i-1)+5   20*(i-1)+9;
        20*(i-1)+9   20*(i-1)+17;
        20*(i-1)+5   20*(i-1)+17;
        20*(i-1)+6   20*(i-1)+11;
        20*(i-1)+11   20*(i-1)+17;
        20*(i-1)+6   20*(i-1)+17;
        20*(i-1)+11   20*(i-1)+22;
        20*(i-1)+17   20*(i-1)+22;
        20*(i-1)+9   20*(i-1)+21;
        20*(i-1)+17   20*(i-1)+21; %40

        20*(i-1)+1    20*(i-1)+5;
        20*(i-1)+2    20*(i-1)+6;
        20*(i-1)+3    20*(i-1)+7;
        20*(i-1)+4    20*(i-1)+8; % add new bars

        ];
end

bar.node_ij_mat=[bar.node_ij_mat;
    20*N+1   20*N+5;
    20*N+2   20*N+6;
    20*N+3   20*N+7;
    20*N+4   20*N+8;
    ];

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);
plots.Plot_Shape_Bar_Number();

%% Define Rotational Spring
for i=1:N    
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        20*(i-1)+5  20*(i-1)+1  20*(i-1)+7  20*(i-1)+2;
        20*(i-1)+1  20*(i-1)+7  20*(i-1)+3  20*(i-1)+8;
        20*(i-1)+7  20*(i-1)+3  20*(i-1)+8  20*(i-1)+4;
        20*(i-1)+3  20*(i-1)+4  20*(i-1)+8  20*(i-1)+6;
        20*(i-1)+2  20*(i-1)+4  20*(i-1)+6  20*(i-1)+8;
        20*(i-1)+4  20*(i-1)+2  20*(i-1)+6  20*(i-1)+1;
        20*(i-1)+2  20*(i-1)+6  20*(i-1)+1  20*(i-1)+5;
        20*(i-1)+6  20*(i-1)+1  20*(i-1)+5  20*(i-1)+7; %8

        20*(i-1)+1  20*(i-1)+7  20*(i-1)+5  20*(i-1)+19;
        20*(i-1)+5  20*(i-1)+7  20*(i-1)+19  20*(i-1)+13;
        20*(i-1)+7  20*(i-1)+13  20*(i-1)+19  20*(i-1)+23;
        20*(i-1)+19  20*(i-1)+23  20*(i-1)+21  20*(i-1)+27;
        20*(i-1)+5  20*(i-1)+10  20*(i-1)+19  20*(i-1)+21;
        20*(i-1)+7  20*(i-1)+5  20*(i-1)+19  20*(i-1)+10;
        20*(i-1)+13  20*(i-1)+19  20*(i-1)+23  20*(i-1)+21;
        20*(i-1)+10  20*(i-1)+19  20*(i-1)+21  20*(i-1)+23; %16

        20*(i-1)+3  20*(i-1)+8  20*(i-1)+7  20*(i-1)+18;
        20*(i-1)+8  20*(i-1)+7  20*(i-1)+18  20*(i-1)+14;
        20*(i-1)+7  20*(i-1)+8  20*(i-1)+18  20*(i-1)+15;
        20*(i-1)+8  20*(i-1)+15  20*(i-1)+18  20*(i-1)+24;
        20*(i-1)+7  20*(i-1)+14  20*(i-1)+18  20*(i-1)+23;
        20*(i-1)+14  20*(i-1)+18  20*(i-1)+23  20*(i-1)+24;
        20*(i-1)+15  20*(i-1)+18  20*(i-1)+24  20*(i-1)+23;
        20*(i-1)+18  20*(i-1)+23  20*(i-1)+24  20*(i-1)+28; %24

        20*(i-1)+4  20*(i-1)+6  20*(i-1)+8  20*(i-1)+20;
        20*(i-1)+6  20*(i-1)+8  20*(i-1)+20  20*(i-1)+16;
        20*(i-1)+8  20*(i-1)+6  20*(i-1)+20  20*(i-1)+12;
        20*(i-1)+8  20*(i-1)+16  20*(i-1)+20  20*(i-1)+24;
        20*(i-1)+6  20*(i-1)+12  20*(i-1)+20  20*(i-1)+22;
        20*(i-1)+24  20*(i-1)+20  20*(i-1)+22  20*(i-1)+12;
        20*(i-1)+16  20*(i-1)+20  20*(i-1)+24  20*(i-1)+22;
        20*(i-1)+26  20*(i-1)+24  20*(i-1)+22  20*(i-1)+20; %32

        20*(i-1)+1  20*(i-1)+5  20*(i-1)+6  20*(i-1)+17;
        20*(i-1)+6  20*(i-1)+5  20*(i-1)+17  20*(i-1)+9;
        20*(i-1)+5  20*(i-1)+6  20*(i-1)+17  20*(i-1)+11;
        20*(i-1)+5  20*(i-1)+9  20*(i-1)+17  20*(i-1)+21;
        20*(i-1)+6  20*(i-1)+17  20*(i-1)+11  20*(i-1)+22;
        20*(i-1)+11  20*(i-1)+17  20*(i-1)+22  20*(i-1)+21;
        20*(i-1)+9  20*(i-1)+17  20*(i-1)+21  20*(i-1)+22;
        20*(i-1)+17  20*(i-1)+21  20*(i-1)+22  20*(i-1)+26; %40
        ];
        
end

rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    20*N+5  20*N+1  20*N+7  20*N+2;
    20*N+1  20*N+7  20*N+3  20*N+8;
    20*N+7  20*N+3  20*N+8  20*N+4;
    20*N+3  20*N+4  20*N+8  20*N+6;
    20*N+2  20*N+4  20*N+6  20*N+8;
    20*N+4  20*N+2  20*N+6  20*N+1;
    20*N+2  20*N+6  20*N+1  20*N+5;
    20*N+6  20*N+1  20*N+5  20*N+7;];


rotNum=size(rot_spr_4N.node_ijkl_mat);
rotNum=rotNum(1);

rot_spr_4N.rot_spr_K_vec=ones(rotNum,1)*100000;

factor=100;
for i=1:N+1
    rot_spr_4N.rot_spr_K_vec((i-1)*40+2)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+2);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+4)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+4);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+6)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+6);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+8)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+8);
end

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

%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

% nodeN=size(node.coordinates_mat,1);
% nr.supp=[(1:nodeN)';zeros(nodeN,1),ones(nodeN,1),zeros(nodeN,1)];

nr.supp=[1 1 1 1;
         2 1 1 1;
         20*N+5 1 1 1;
         20*N+6 1 1 1;];

force=50000; 
step=100;
Uhis=[];

A=bar.A_vec(:); % m^2
E=bar.E_vec(:); % Pa
K =1.0; % Assume Pin-Pin for every bar               
Lc=K*bar.L0_vec;

r=0.0287*ones(size(A)); %0.5*sqrt(bar.A_vec./pi);
Fy=3.45e8;  % Pa

nb=size(A,1);

for k=1:step

    nr.load=[N*10+1,0,0,-force*k/4;
             N*10+2,0,0,-force*k/4;
             N*10+5,0,0,-force*k/4;
             N*10+6,0,0,-force*k/4];
    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=1*10^-5;

    Uhis(k,:,:)=squeeze(nr.Solve());

    truss_strain=bar.Solve_Strain(node,squeeze(Uhis(k,:,:)));
    Nforce=(truss_strain).*(bar.E_vec).*(bar.A_vec);

    [passVec,modeStrList]=CheckAllBarDesignAISC(Nforce,A,E,Lc,r,Fy);

    if sum(passVec)==nb
    else
        
        fprintf('Failure happended \n');

        % Plot the deformed shape when failure happen
        plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)))
        
        % Find the internal force when failure occur
        truss_strain=bar.Solve_Strain(node,squeeze(Uhis(end,:,:)));
        internal_force=(truss_strain).*(bar.E_vec).*(bar.A_vec);
        
        % Plot failure stress
        bar_stress=(truss_strain).*(bar.E_vec);
        plots.Plot_Shape_Bar_Stress(bar_stress)
        
        % Find total bar length
        barLtotal=sum(bar.L0_vec);
        
        % Find Stiffness
        Uaverage=-mean(squeeze(Uhis(end,[N*10+1,N*10+2,N*10+5,N*10+6],3)));
        Kstiff=step*force/Uaverage;
        
        % Find the relationship betweeen the bar internal forces and load
        loadatfail=force*step;
        fprintf('Failure load is %d kN \n',  loadatfail/1000);
        fprintf('Total bar length is %d m \n',  barLtotal);
        fprintf('Stiffness is %d N/m \n',  Kstiff);
        fprintf('Total bars: %d\n', barNum);

        plots.Plot_Shape_Bar_Failure(passVec)

        break
    end

    if k==step
        fprintf('No Failure Happen \n', barNum);
    end
end

% % Plot the deformed shape when failure happen
% plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)))
% 
% % Find the internal force when failure occur
% truss_strain=bar.Solve_Strain(node,squeeze(Uhis(end,:,:)));
% internal_force=(truss_strain).*(bar.E_vec).*(bar.A_vec);
% 
% % Plot failure stress
% bar_stress=(truss_strain).*(bar.E_vec);
% plots.Plot_Shape_Bar_Stress(bar_stress)
% 
% % Find total bar length
% barLtotal=sum(bar.L0_vec);
% 
% % Find Stiffness
% Uaverage=-mean(squeeze(Uhis(end,[N*10+1,N*10+2,N*10+5,N*10+6],3)));
% Kstiff=step*force/Uaverage;
% 
% % Find the relationship betweeen the bar internal forces and load
% loadatfail=force*step;
% fprintf('Failure load is %d kN \n',  loadatfail/1000);
% fprintf('Total bar length is %d m \n',  barLtotal);
% fprintf('Stiffness is %d N/m \n',  Kstiff);
% fprintf('Total bars: %d\n', barNum);
% 
% loadEff=loadatfail/W_bar


% %% Evaluate Member
% AxialForce=internal_force(:); % N
% A=bar.A_vec(:); % m^2
% E=bar.E_vec(:); % Pa
% nb=numel(A);
% 
% % 1) effective length KL
% if ~isfield(bar,'L0_vec')||isempty(bar.L0_vec)
%     L0_vec=zeros(nb,1);
%     for k=1:nb
%         n1=bar.node_ij_mat(k,1);
%         n2=bar.node_ij_mat(k,2);
%         L0_vec(k)=norm(assembly.node.coordinates_mat(n1,:) - ...
%                          assembly.node.coordinates_mat(n2,:));
%     end
% else
%     L0_vec=bar.L0_vec(:);
% end
% K =1.0;                
% Lc=K.*L0_vec;
% 
% % 2) r, r = 0.5*sqrt(A/pi)
% r=0.5*sqrt(A./pi);
% r=max(r,1e-9); % prevent division by zero
% 
% % 3) yield stress
% Fy=250e6;             % Pa（need to be changed）
% 
% % 4) evaluate
% passYN=false(nb,1);
% util=NaN(nb,1);
% modeStr=cell(nb,1);
% Pn=NaN(nb,1);
% slender=NaN(nb,1);   % KL/r
% Fe=NaN(nb,1);   % Euler stress (Pa)
% Fcr=NaN(nb,1);   % Critical stress per AISC (Pa)
% 
% tiny=1e-12;
% for i=1:nb
%     Ni=AxialForce(i);
%     Ai=A(i);
%     Ei=E(i);
%     Lci=Lc(i);
%     ri=r(i);
% 
%     if Ni>0
%         % in tension 
%         Pn_i=Fy*Ai;
%         modeStr{i}='Tension-Yield';
%         slender(i)=NaN;  
%         Fe(i)=NaN; 
%         Fcr(i)=NaN;   % don't calculate buckling in tension
%     else
%         % in compression
%         slender=Lci/ri;
%         Fe=(pi^2*Ei)/(slender^2);
%         lambda_lim=4.71*sqrt(Ei/Fy);
%         if slender<=lambda_lim
%             Fcr=(0.658)^(Fy/Fe)*Fy;   % inelastic buckling
%         else
%             Fcr=0.877*Fe;             % elastic bukling
%         end
%         Pn_i=Fcr*Ai;
%         modeStr{i}='Compression-Buckling';
%     end
% 
%     Pn(i)=Pn_i;
%     util(i)=abs(Ni)/max(Pn_i,tiny);
%     passYN(i)=util(i)<=1.0;
% end

% % 5) print
% fprintf('Bars passed: %d / %d\n', sum(passYN), nb);
% 
% [util_sorted, idx]=sort(util, 'descend');
% topk=min(10, nb);
% fprintf('Worst %d bars (by utilization):\n', topk);
% for ii=1:topk
%     b=idx(ii);
%     fprintf('#%d: N=%.2f kN, util=%.3f, mode=%s, Pn=%.2f kN\n', ...
%         b, AxialForce(b)/1e3, util_sorted(ii), modeStr{b}, Pn(b)/1e3);
% end


% function [pass,modeStr]=CheckBarDesignAISC(Ni,Ai,Ei,Lci,ri,Fy)
% 
%     tiny=1e-12;
% 
%     if Ni>0
%         % in tension 
%         Pn_i=Fy*Ai;
%         modeStr='Tension-Yield';
%         slender=NaN;  
%         Fe=NaN; 
%         Fcr=NaN;   % don't calculate buckling in tension
%     else
%         % in compression
%         slender=Lci/ri;
%         Fe=(pi^2*Ei)/(slender^2);
%         lambda_lim=4.71*sqrt(Ei/Fy); % threshold for elastic/inelastic buckling
% 
%         if slender<=lambda_lim
%             Fcr=(0.658)^(Fy/Fe)*Fy;   % inelastic buckling
%         else
%             Fcr=0.877*Fe;             % elastic bukling
%         end
%         Pn_i=Fcr*Ai;
%         modeStr='Compression-Buckling';
%     end
% 
%     Pn=Pn_i;
%     util=abs(Ni)/max(Pn_i,tiny);
%     pass=util<=1.0;
% 
% end
% 
% function [passVec,modeStrList]=CheckAllBarDesignAISC(Nforce,A,E,Lc,r,Fy)
%     nbar=size(A,1);
%     passVec=zeros(nbar,1);
%     modeStrList={};
%     for i=1:nbar
% 
%         [pass,modeStr]=CheckBarDesignAISC(Nforce(i),A(i),E(i),Lc(i),r(i),Fy);
%         passVec(i)=pass;
%         modeStrList{i}=modeStr;
%     end
% end