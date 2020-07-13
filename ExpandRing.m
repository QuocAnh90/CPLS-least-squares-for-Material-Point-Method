function [results] = ExpandRing(Version, resolution)

% clear all 
% close all

addpath('SubFunctions');
% Unit
% Newton - seconds - metre

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'MPM'
% CPDI: 'CPDI'

% Please remember to change the output video file!!
% Version = 'IMLS_CPDI';

%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
CModel = 'Neo_Hookean_Elastic';

%% Material porperties
E                       = 10000000          ;               % Young modulus of solid
psp                     = 1000.0            ;               % solid density
nu                      = 0.0           ;                   % Poison ratio
g                       = 0.0            ;                  % gravity acceleration
Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);
CModel_parameter = [E,nu];

%% Structured Grid input
% resolution              = 2;
le(1)                   = 0.1/resolution;                    % size of element in X direction
le(2)                   = 0.1/resolution;                    % size of element in Y direction
NN(1)                   = 12*resolution+1;                       % number of nodes in X direction
NN(2)                   = 12*resolution+1;                       % number of nodes in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

%% Time
C                       = sqrt(E/psp);
ftime                   = 0.02;
dt                      = 0.4*le(1)/C;
ndt                     = round(ftime/dt) +1;
t                       = 0;
n                       = 1;

spElems_boundary = 0;
p_sc = zeros(cellCount,1);
color_node = zeros(nodeCount,1);
gridmomentum=[];
particlemomentum=[];

%% Boundary nodes
% Boundary coordination
x_min = 2*le(1);
x_max = (NN(1) - 3)*le(1);
y_min = 2*le(2);
y_max = (NN(2) - 3)*le(2);
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);

%% Particle generation
particle_per_cell       = 4;
spCount                 = 10*10*resolution*resolution*particle_per_cell;
lp(1)                   = le(1)/sqrt(particle_per_cell);                                 % size of particle in X direction
lp(2)                   = le(2)/sqrt(particle_per_cell);                                % size of particle in Y direction


nPar1                     = 16* resolution;
nPar2                     = 2 * 2 * resolution;
spCount                 = nPar1*nPar2;
x_sp                    = zeros(spCount,2);
d_sp                    = zeros(spCount,2);
V_sp                    = zeros(spCount,1);
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);

color_node              = zeros(nodeCount,1);

% Ring distribution
for i = 1:nPar1
    ttemp=(i-0.5)*0.5*pi/nPar1;
    for j = 1:nPar2
        rtemp=0.4+(j-0.5)*0.2/nPar2;
        x_sp((i-1)*nPar2+j,:)=[2*le(1)+rtemp*cos(ttemp);2*le(2)+rtemp*sin(ttemp)];
        V_sp((i-1)*nPar2+j)=0.2/nPar2*0.5*pi/nPar1*rtemp;
        r1_sp((i-1)*nPar2+j,:)=0.25*pi/nPar1*rtemp*[sin(ttemp);-cos(ttemp)];
        r2_sp((i-1)*nPar2+j,:)=0.1/nPar2*[cos(ttemp);sin(ttemp)];
    end
end

% %% Square distributions
% sp=1;
% d_sp                    = zeros(spCount,2);
% V_sp                    = zeros(spCount,1);
% r1_sp                   = zeros(spCount,2);
% r2_sp                   = zeros(spCount,2);
% while sp<spCount+0.0001
%     for i=1:10*resolution*1*sqrt(particle_per_cell)
%         for j=1:10*resolution*sqrt(particle_per_cell)
%             x_sp(sp,1:2)= [1*le(1)+0.5*lp(1,1)+(j-1)*lp(1,1) 1*le(2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
%             sp=sp+1;
%         end
%     end  
% end
% sp=1;
% while sp<spCount+1
%         x_center = x_sp(sp,1)-le(1);
%         y_center = x_sp(sp,2)-le(2);
%         Radius = 0.6;
%         
%      if (x_center^2+y_center^2)>Radius^2 
%     x_sp(sp,:)=[];
%     r1_sp(sp,:)=[];
%     r2_sp(sp,:)=[];
%     V_sp(sp)=[];
%     spCount=spCount-1;
%     sp=sp-1;
%      end
%     sp=sp+1;
% end
% 
% sp=1;
% while sp<spCount+1
%         x_center = x_sp(sp,1)-le(1);
%         y_center = x_sp(sp,2)-le(2);
%         Radius = 0.4;
%         
%      if (x_center^2+y_center^2)<Radius^2 
%     x_sp(sp,:)=[];
%     r1_sp(sp,:)=[];
%     r2_sp(sp,:)=[];
%     V_sp(sp)=[];
%     spCount=spCount-1;
%     sp=sp-1;
%      end
%     sp=sp+1;
% end
% 
% % Initial volume
% for sp = 1:spCount
%     r1_sp(sp,:) = [lp(1,1)/2 0];
%     r2_sp(sp,:) = [0 lp(1,2)/2];
%     V_sp(sp)    = 4*abs(r1_sp(sp,1)*r2_sp(sp,2)-r1_sp(sp,2)*r2_sp(sp,1)); 
% end

V_spo = V_sp;
r10_sp = r1_sp;
r20_sp = r2_sp;

x_corner1 = x_sp - r1_sp - r2_sp;
x_corner2 = x_sp + r1_sp - r2_sp;
x_corner3 = x_sp + r1_sp + r2_sp;
x_corner4 = x_sp - r1_sp + r2_sp;

x_data = cell(spCount,5);
% boundary_corner = [];

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
% stress_data = cell(spCount,5);
% density_data = cell(spCount,5);
% body_data = cell(spCount,5);

for sp = 1:spCount
    velocity_data{sp,1} = body_ringexpansion(x_sp(sp,1)-2*le(1),x_sp(sp,2)-2*le(2),0,2);
    velocity_data{sp,2} = body_ringexpansion(x_corner1(sp,1)-2*le(1),x_corner1(sp,2)-2*le(2),0,2);    
    velocity_data{sp,3} = body_ringexpansion(x_corner2(sp,1)-2*le(1),x_corner2(sp,2)-2*le(2),0,2);
    velocity_data{sp,4} = body_ringexpansion(x_corner3(sp,1)-2*le(1),x_corner3(sp,2)-2*le(2),0,2);
    velocity_data{sp,5} = body_ringexpansion(x_corner4(sp,1)-2*le(1),x_corner4(sp,2)-2*le(2),0,2);
end

%% Particle variables
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
d_sp                    = zeros(spCount,2);                     % displacement
b_sp                    = zeros(spCount,2);                     % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
v_ssp                   = zeros(spCount,2);                     % velocty
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation

%% Initial condition
% Gradient deformation
for sp = 1:spCount
    F_sp{sp} = [1 0; 0 1];
end
m_sp                    = psp * V_sp;                           % mass

% Generate the velocity
for sp=1:spCount
    xptemp=[x_spo(sp,1)-2*le(1);x_spo(sp,2)-2*le(2)];
    v_ssp(sp,:) = body_ringexpansion(xptemp(1),xptemp(2),0,2);
end

% % start the algorithm
% % video
% timestep = 20;     % number of frame to save
% r=timestep/10;      % number of frame per second video ~200s
% 
% writerObj2           = VideoWriter('RingExpand.avi');
% writerObj2.FrameRate = r;    % number of frame per second
% open(writerObj2);
% 
%     for tt = 1:timestep
%     ft              = ftime/timestep*tt;
    ft=ftime;
 while t<ft+0.000000000001
% for n=1:400
     Version
     resolution
     t

     % Methods of Manufactued solution
     for sp=1:spCount      
     % Body force applied to the particles
        xptemp=[x_spo(sp,1)-2*le(1);x_spo(sp,2)-2*le(2)];
        b_sp(sp,:) = body_ringexpansion(xptemp(1),xptemp(2),(n-1)*dt,1);
     end
     
     switch Version
         case 'MPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MPM_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);
   
         case 'CPDI'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp] = CPDI_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
    
        case 'IMLS_CPDI'   
        [gridmomentum,particlemomentum,k,~,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,...
        r1_sp,r2_sp,L_sp,velocity_data,x_data]...
        = IMLS_CPDI_solver_test(gridmomentum,particlemomentum,...
        velocity_data,x_data,...
        CModel,CModel_parameter,...
        nodeCount,spCount,cellCount...
        ,x_sp,x_spo,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
    
     end
    
     % Methods of Manufactued solution
     % Compute the analytical solutions
     for sp=1:spCount
        xptemp          = [x_spo(sp,1)-2*le(1);x_spo(sp,2)-2*le(2)];
        % Displacement
        dis_ana(sp,:)   = body_ringexpansion(xptemp(1),xptemp(2),n*dt,3);        
        % Velocity
        vel_ana(sp,:)   = body_ringexpansion(xptemp(1),xptemp(2),n*dt,2);        
        % Stress
        stress_ana{sp} = body_ringexpansion(xptemp(1),xptemp(2),n*dt,5);        
     end
%          StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,0,v_ssp,spCount,r1_sp,r2_sp);

%       Update time and step 
        t = t+dt;
        n = n+1;
 end

% %  Plot the result
%     StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,color_node,v_ssp,spCount,spElems_boundary,cellCount,nodeCount,r1_sp,r2_sp,p_sc,LOCC);  
%     frame2 = getframe(StressProfile1);
%     writeVideo(writerObj2,frame2);
%     end
%     close(writerObj2);
    
     %% Compute error
    x=0;    v=0;    s=0;
    deviation_X = zeros(spCount,2);    
    deviation_V = zeros(spCount,2);    
    deviation_S = cell(spCount,1);
    
    % Loop all particles to store the RMS
    for sp=1:spCount      
    deviation_X(sp,:) = dis_ana(sp,:)-d_sp(sp,:);
    x=x + norm(deviation_X(sp,:))^2;
    
    deviation_V(sp,:) = vel_ana(sp,:)-v_ssp(sp,:);
    v=v + norm(deviation_V(sp,:))^2;
    
    deviation_S{sp} = stress_ana{sp}-[s_sp(sp,1) s_sp(sp,3);s_sp(sp,3) s_sp(sp,2)];
    s=s + sqrt(deviation_S{sp}(1,1)^2 + deviation_S{sp}(2,1)^2 + deviation_S{sp}(2,2)^2);
    end    
    perror_X = sqrt(x/spCount); perror_V = sqrt(v/spCount); perror_S = sqrt(s/spCount);   
    results = [perror_X perror_V perror_S];