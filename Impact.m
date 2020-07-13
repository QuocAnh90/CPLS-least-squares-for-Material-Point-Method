% function results = Impact(Version,resolution)

addpath('SubFunctions');
clear all;
close all;

% Unit
% Newton - seconds - metre

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
Version = 'IMLS_CPDI';

%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
CModel = 'Neo_Hookean_Elastic';

%% Material porperties
E                       = 10000000          ;                % Young modulus of solid
psp                     = 1000.0            ;               % solid density
nu                      = 0.0           ;                   % Poison ratio
g                       = 0.0            ;                  % gravity acceleration
Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);
CModel_parameter = [E,nu];

%% Structured Grid input
resolution              = 200;
NN(1)                   = resolution + 5;                  % number of nodes in X direction
NN(2)                   = 6;                  % number of nodes in Y direction
le(1)                   = 1/resolution;                    % size of element in X direction
le(2)                   = 1/resolution;                    % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

%% Time
C                       = sqrt(E/psp);
ftime                   = 0.75*1/C;
dt                      =  0.001*1/C;
ndt                     = round(ftime/dt) +1;
t                       = 0;

%% Boundary nodes
% Boundary coordination
x_min = 2*le(1);
x_max = 10;
y_min = 2*le(2);
y_max = 3*le(2);
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

%% Particle generation
particle_per_cell       = 4;
spCount                 = 1*resolution*particle_per_cell;
lp(1)                   = le(1)/sqrt(particle_per_cell);                                 % size of particle in X direction
lp(2)                   = le(2)/sqrt(particle_per_cell);                                % size of particle in Y direction
x_sp                    = zeros(spCount,2);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:1*1*sqrt(particle_per_cell)
        for j=1:resolution*sqrt(particle_per_cell)
            x_sp(sp,1:2)= [2*le(2)+0.5*lp(1,1)+(j-1)*lp(1,1) 2*le(2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
            sp=sp+1;
        end
    end  
end
% x_sp: Vector, position of MPs
% spCount: total number of MPs

%% Plot initial condition
initial_figure = Plot_Initial(x_sp,LOC,le);

%% Particle variables
dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
d_sp                    = zeros(spCount,2);                     % displacement
b_sp                    = zeros(spCount,2);                     % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
a_ssp                   = zeros(spCount,2);                     % acceleration
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);

% Plot
spElems_boundary        = 0;
p_sc                    = zeros(cellCount,1);
color_node              = zeros(nodeCount,1);
b_sc                    = zeros(cellCount,2);

%% Initial condition
% Gradient deformation
for sp = 1:spCount
    r1_sp(sp,:) = [lp(1,1)/2 0];
    r2_sp(sp,:) = [0 lp(1,2)/2];
    F_sp{sp} = [1 0; 0 1];
end
r10_sp = r1_sp;
r20_sp = r2_sp;
V_sp                    = zeros(spCount,1);
for sp=1:spCount
V_sp(sp)                = 4*abs(r1_sp(sp,1)*r2_sp(sp,2)-r1_sp(sp,2)*r2_sp(sp,1)); 
end
V_spo                   = V_sp;
m_sp                    = psp * V_sp;                           % mass

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
    velocity_data{sp,1} = [0;0];
    velocity_data{sp,2} = [0;0];
    velocity_data{sp,3} = [0;0];
    velocity_data{sp,4} = [0;0];
    velocity_data{sp,5} = [0;0];
    
    acceleration_data{sp,1} = [0;0];
    acceleration_data{sp,2} = [0;0];
    acceleration_data{sp,3} = [0;0];
    acceleration_data{sp,4} = [0;0];
    acceleration_data{sp,5} = [0;0];
end

% start the algorithm
% video
timestep = 100;     % number of frame to save
r=timestep/10;      % number of frame per second video ~200s

writerObj2           = VideoWriter('aligned1Dx.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:timestep
    ft              = ftime/timestep*tt;
%     ft=ftime;
 while t<ft+0.000001      
     t
     
     ptraction_sp(:)=0;
     % Traction
if t < 0.5*(1)/C && t>0*(1)/C
    for p = 1:spCount
    if x_spo(p,1) > 1+le(1)
        ptraction_sp(p,1) = 1;
    end
    end
end

     switch Version
         case 'MPM'
%         [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MPM_solver(CModel,CModel_parameter,...
%         nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
%         nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);
        
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,a_ssp] = MPM_solver_alpha(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,lp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,a_ssp);

         case 'CPDI'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp] = CPDI_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount...
        ,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
    
%         [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp,a_ssp]...
%         = CPDI_solver_alpha(CModel,CModel_parameter,...
%         nodeCount,spCount,cellCount,...
%         x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,...
%         v_ssp,s_sp,m_sp,p_sp,nfbcx,nfbcy,fbcx,fbcy,...
%         F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp,a_ssp);
    
        case 'IMLS_CPDI'
%         [DensityGradient,StressGradient,p_sc,spElems_boundary,CellBoundary_node,color_node,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,...
%         r1_sp,r2_sp,L_sp,body_data,velocity_data,stress_data,density_data,x_data]...
%         = IMLS_CPDI_solver(boundary_corner,...
%         body_data,velocity_data,stress_data,density_data,x_data,...
%         CModel,CModel_parameter,...
%         nodeCount,spCount,cellCount...
%         ,x_sp,x_spo,le,NN,LOC,LOCC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
%         nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
    
%         [color_node,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,...
%         r1_sp,r2_sp,L_sp,velocity_data,x_data]...
%         = IMLS_CPDI_solver_test(...
%         velocity_data,x_data,...
%         CModel,CModel_parameter,...
%         nodeCount,spCount,cellCount...
%         ,x_sp,x_spo,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
%         nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
% 
        [color_node,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,...
        r1_sp,r2_sp,L_sp,A_sp,velocity_data,acceleration_data,x_data]...
        = IMLS_CPDI_solver_test1(...
        velocity_data,acceleration_data,x_data,...
        CModel,CModel_parameter,...
        nodeCount,spCount,cellCount...
        ,x_sp,x_spo,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,a_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
    
     end
     
 % Update time and step 
 t = t+dt;
 end


 %% Plot the result
% StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,color_node,v_ssp,spCount,spElems_boundary,cellCount,nodeCount,r1_sp,r2_sp,p_sc,LOCC);  
    
% for p=1:spCount
%     x_sp(p,1) = x_sp(p,1)-2*le(1);
% end

x1=[0 0.25];y1=[0 0];x2=[0.25 0.25];y2=[0 1];x3=[0.25 0.75];y3=[1 1];x4=[0.75 0.75];y4=[1 0];x5=[0.75 1];y5=[0 0];
    StressProfile1=figure;
    plot(x_sp(1:spCount/2,1),s_sp(1:spCount/2,1),'.',x1,y1,'g',x2,y2,'g',x3,y3,'g',x4,y4,'g',x5,y5,'g');
    line(x_sp(1:spCount/2,1),s_sp(1:spCount/2,1));
%     axis([0 1 -0.4 1.4])
    ylabel('Stress (Pa)'); % label for y axis
    xlabel('Length (s)'); % label for x axis
    legend('CPDI','analytical solution','Location','southwest');
    set(gca,'fontsize', 14)

% for p=1:spCount
%     x_sp(p,1) = x_sp(p,1)+2*le(1);
% end

frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    close(writerObj2);
