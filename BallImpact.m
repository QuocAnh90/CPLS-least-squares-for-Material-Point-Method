% close all;
clear all;
tic;
addpath('SubFunctions/');

% Unit
% Newton - seconds - metre

%% Please select the versions of MPM
% Original MPM: 'MPM'
% CPDI: 'CPDI'
% Please remember to change the output video file!!
Version = 'CPDI';

%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
CModel = 'Linear_Elastic';

%% Material porperties
E                       = 1000           ;                  % Young modulus of solid
psp                     = 1000.0         ;                  % solid density
nu                      = 0.3            ;                  % Poison ratio
g                       = 0.0            ;                  % gravity acceleration

CModel_parameter = [E,nu];

%% Structured Grid input
NN(1)                 = 25;                               % number of nodes in X direction
NN(2)                 = 25;                               % number of nodes in Y direction
le(1)                 = 0.05;                              % size of element in X direction
le(2)                 = le(1,1);                          % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

%% Time
ftime                   = 3.5;                               % Final time
dt                      = 0.005;                            % Time step
ndt                     = round(ftime/dt) +1;
t                       = 0;

% Plot
spElems_boundary = 0;
p_sc = zeros(cellCount,1);
color_node = zeros(nodeCount,1);
gridmomentum=[];
particlemomentum=[];
time=[];
NNN =[];

%% Boundary nodes
% Boundary coordination
x_min = -100;
x_max = 100;
y_min = -100;
y_max = 100;
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

%% Particle generation
[spCount,x_sp,lp]=Two_Particle_Ball(le,0.3,0.3,0.8,0.8,0.2);
% x_sp: Vector, position of MPs
% spCount: total number of MPs

%% Plot initial condition
initial_figure = Plot_Initial(x_sp,LOC,le);

%% Particle variables
dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
d_sp                    = zeros(spCount,2);                     % displacement
V_sp                    = dparticle * ones(spCount,1);          % volumn
V_spo                   = V_sp;                                 % initial volumn
m_sp                    = psp * V_sp;                           % mass
b_sp                    = zeros(spCount,2);                     % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);
a_ssp                   = zeros(spCount,2);

%% Initial condition
% Gradient deformation
for sp = 1:spCount
    r1_sp(sp,:) = [lp(1,1)/2 0];
    r2_sp(sp,:) = [0 lp(1,2)/2];
    F_sp{sp} = [1 0; 0 1];
end
r10_sp = r1_sp;
r20_sp = r2_sp;

% Velocity condition
for sp=1:spCount
    if x_sp(sp,1)<0.55
    v_ssp(sp,:) = [sqrt(2)/20 sqrt(2)/20];
    elseif x_sp(sp,1)>0.55
    v_ssp(sp,:) = [-sqrt(2)/20 -sqrt(2)/20];
    end
end

% Generate for IMLS_CPDI
x_corner1 = x_sp - r1_sp - r2_sp;
x_corner2 = x_sp + r1_sp - r2_sp;
x_corner3 = x_sp + r1_sp + r2_sp;
x_corner4 = x_sp - r1_sp + r2_sp;

x_data = cell(spCount,5);
for sp = 1:spCount
    x_data{sp,1} = x_sp(sp,:)';
    x_data{sp,2} = x_corner1(sp,:)';
    x_data{sp,3} = x_corner2(sp,:)';
    x_data{sp,4} = x_corner3(sp,:)';
    x_data{sp,5} = x_corner4(sp,:)';
end

x_data_initial = x_data;

%% Build velocity data for IMLS
velocity_data = cell(spCount,5);
acceleration_data = cell(spCount,5);

for sp = 1:spCount  
    if x_sp(sp,1)<0.55
    velocity_data{sp,1} = [sqrt(2)/20;sqrt(2)/20];
    velocity_data{sp,2} = [sqrt(2)/20;sqrt(2)/20];
    velocity_data{sp,3} = [sqrt(2)/20;sqrt(2)/20];
    velocity_data{sp,4} = [sqrt(2)/20;sqrt(2)/20];
    velocity_data{sp,5} = [sqrt(2)/20;sqrt(2)/20];
    elseif x_sp(sp,1)>0.55
    velocity_data{sp,1} = [-sqrt(2)/20;-sqrt(2)/20];
    velocity_data{sp,2} = [-sqrt(2)/20;-sqrt(2)/20];
    velocity_data{sp,3} = [-sqrt(2)/20;-sqrt(2)/20];
    velocity_data{sp,4} = [-sqrt(2)/20;-sqrt(2)/20];
    velocity_data{sp,5} = [-sqrt(2)/20;-sqrt(2)/20];
    end
    
    acceleration_data{sp,1} = [0;0];
    acceleration_data{sp,2} = [0;0];
    acceleration_data{sp,3} = [0;0];
    acceleration_data{sp,4} = [0;0];
    acceleration_data{sp,5} = [0;0];
end


% Data
aM = [];
M = [];
time =[];
kinetic=[];
strain=[];
total=[];

%% start the algorithm
% video
timestep = 100;     % number of frame to save
r=timestep/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('BallImpact.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:timestep
    ft              = ftime/timestep*tt;
%     ft=ftime;
 while t<ft+0.0000000001      
     t
     
    switch Version
         case 'MPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MPM_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);

         case 'CPDI'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp] = CPDI_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount...
        ,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
    
        case 'IMLS_CPDI'
        [gridmomentum,particlemomentum,k,color_node,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,...
        r1_sp,r2_sp,L_sp,velocity_data,x_data,nvelo_si]...
        = IMLS_CPDI_solver_test(gridmomentum,particlemomentum,...
        velocity_data,x_data,...
        CModel,CModel_parameter,...
        nodeCount,spCount,cellCount...
        ,x_sp,x_spo,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
    end
 
time = [time;t];
     
    % Energy
    KineticEnergy = 0;
    StrainEnergy = 0;
    Total = 0;
    
    for p = 1:spCount
        KineticEnergy = KineticEnergy + m_sp(p)*(v_ssp(p,1)^2+v_ssp(p,2)^2)/2;
        StrainEnergy = StrainEnergy + 0.5 * V_sp(p)*(s_sp(p,:)*s_sp(p,:)')/E;
    end
    
    Total = KineticEnergy + StrainEnergy;
    kinetic=[kinetic;KineticEnergy];
    strain=[strain;StrainEnergy];
    total=[total;Total];    
     
 % Update time and step 
 t = t+dt;
 end
% 
%  %% Plot the result
    StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,color_node,v_ssp,spCount,spElems_boundary,cellCount,nodeCount,r1_sp,r2_sp,p_sc,LOCC);  
% 
    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    

    figure
    plot(time,kinetic,time,strain,time,total)
    legend('kinetic','strain','total','Location','south');
    
    