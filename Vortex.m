function [results] = Vortex(Version, resolution)

% close all;
% clear all;
% % tic;

addpath('SubFunctions');

% Unit
% Newton - seconds - metre

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'MPM'
% CPDI: 'CPDI'
% CPDI: 'IMLS_CPDI'

% Please remember to change the output video file!!
% Version = 'CPDI';
A = 1;
%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
CModel = 'Neo_Hookean_Elastic';

%% Material porperties
E                       = 1000          ;                  % Young modulus of solid
psp                     = 1000.0        ;                  % solid density
nu                      = 0.3           ;                  % Poison ratio
g                       = 0.0           ;                  % gravity acceleration
Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);
CModel_parameter = [E,nu];

%% Structured Grid input
% resolution              = 2;
NN(1)                   = 16*resolution+1;                               % number of nodes in X direction
NN(2)                   = 16*resolution+1;                               % number of nodes in Y direction
le(1)                   = 0.25/resolution;                              % size of element in X direction
le(2)                   = 0.25/resolution;                              % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction
center = max(LOC(:,1))/2;

%% Time
C                       = sqrt(E/psp);
ftime                   = 1.0;
dt                      = 0.04*le(1)/C;
ndt                     = round(ftime/dt) +1;
t                       = 0;
n                       = 1;

% Plot
spElems_boundary = 0;
p_sc = zeros(cellCount,1);
color_node = zeros(nodeCount,1);
gridmomentum=[];
particlemomentum=[];

% Particle generation
% % Square distrobution
% particle_per_cell       = 4;
% spCount                 = 14*resolution*14*resolution*particle_per_cell;
% lp(1)                   = le(1)/sqrt(particle_per_cell);                                 % size of particle in X direction
% lp(2)                   = le(2)/sqrt(particle_per_cell);                                % size of particle in Y direction
% x_sp                    = zeros(spCount,2);
% V_sp                    = zeros(spCount,1);
% r1_sp                   = zeros(spCount,2);
% r2_sp                   = zeros(spCount,2);
% 
% sp=1;
% while sp<spCount+0.0001
%     for i=1:14*resolution*sqrt(particle_per_cell)
%         for j=1:14*resolution*sqrt(particle_per_cell)
%             x_sp(sp,1:2)= [2*le(2)+0.5*lp(1,1)+(j-1)*lp(1,1) 2*le(2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
%             sp=sp+1;
%         end
%     end  
% end
% 
% sp=1;
% while sp<spCount+1
%         x_center = x_sp(sp,1)-center;
%         y_center = x_sp(sp,2)-center;
%         Radius = 1.25;
%         
%      if (x_center^2+y_center^2)>Radius^2 
%     x_sp(sp,:)=[];
%     V_sp(sp,:)=[];
%     r1_sp(sp,:)=[];
%     r2_sp(sp,:)=[];
%     spCount=spCount-1;
%     sp=sp-1;
%      end
%     sp=sp+1;
% end
% 
% sp=1;
% while sp<spCount+1
%         x_center = x_sp(sp,1)-center;
%         y_center = x_sp(sp,2)-center;
%         Radius = 0.75;
%         
%      if (x_center^2+y_center^2)<Radius^2 
%     x_sp(sp,:)=[];
%     V_sp(sp,:)=[];
%     r1_sp(sp,:)=[];
%     r2_sp(sp,:)=[];
%     spCount=spCount-1;
%     sp=sp-1;
%      end
%     sp=sp+1;
% end
% 
% for sp = 1:spCount
%     r1_sp(sp,:) = [lp(1,1)/2 0];
%     r2_sp(sp,:) = [0 lp(1,2)/2];
%     V_sp(sp)    = 4*abs(r1_sp(sp,1)*r2_sp(sp,2)-r1_sp(sp,2)*r2_sp(sp,1));
% end

% Ring distribution
nPar1                   = 45* resolution;
nPar2                   = 4 * resolution;
spCount                 = nPar1*nPar2;
x_sp                    = zeros(spCount,2);
V_sp                    = zeros(spCount,1);
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);

% Ring distribution
for i = 1:nPar1
    ttemp=(i-0.5)*2*pi/nPar1;
    for j = 1:nPar2
        rtemp=0.75+(j-0.5)*0.5/nPar2;
        x_sp((i-1)*nPar2+j,:)=[center+rtemp*cos(ttemp);center+rtemp*sin(ttemp)];
        V_sp((i-1)*nPar2+j)=0.5/nPar2*2*pi/nPar1*rtemp;
        r1_sp((i-1)*nPar2+j,:)=1*pi/nPar1*rtemp*[sin(ttemp);-cos(ttemp)];
        r2_sp((i-1)*nPar2+j,:)=0.25/nPar2*[cos(ttemp);sin(ttemp)];
    end
end

V_spo = V_sp;
r10_sp = r1_sp;
r20_sp = r2_sp;

% Corner position
x_corner1 = x_sp - r1_sp - r2_sp;
x_corner2 = x_sp + r1_sp - r2_sp;
x_corner3 = x_sp + r1_sp + r2_sp;
x_corner4 = x_sp - r1_sp + r2_sp;

% Compute x_data and boundary corners
x_data = cell(spCount,5);
boundary_corner = [];

for sp = 1:spCount
    x_data{sp,1} = x_sp(sp,:)';
    x_data{sp,2} = x_corner1(sp,:)';
    x_data{sp,3} = x_corner2(sp,:)';
    x_data{sp,4} = x_corner3(sp,:)';
    x_data{sp,5} = x_corner4(sp,:)';
end

x_data_initial = x_data;


%% Boundary nodes
% Boundary coordination
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

nfbcx                   = 0              ;  % initial number of fixed nodes in x direction
nfbcy                   = 0              ;  % initial number of fixed nodes in y direction
fbcx = []; fbcy = [];                       % vector store the index of boundary nodes

%% Plot initial condition
initial_figure = Plot_Initial(x_sp,LOC,le);
hold on 

for sp = 1:spCount
%     Qxx{sp} = [x_sp(sp,1) x_sp(sp,1)+r1_sp(sp,1)];
%     Qxy{sp} = [x_sp(sp,2) x_sp(sp,2)+r1_sp(sp,2)];
%     Qyx{sp} = [x_sp(sp,1) x_sp(sp,1)+r2_sp(sp,1)];
%     Qyy{sp} = [x_sp(sp,2) x_sp(sp,2)+r2_sp(sp,2)];
%     plot (Qxx{sp},Qxy{sp},'r')
%     hold on 
%     plot (Qyx{sp},Qyy{sp},'r')
%     hold on 
 
 recx{sp} = [x_corner1(sp,1) x_corner2(sp,1) x_corner3(sp,1) x_corner4(sp,1) x_corner1(sp,1)];
 recy{sp} = [x_corner1(sp,2) x_corner2(sp,2) x_corner3(sp,2) x_corner4(sp,2) x_corner1(sp,2)];
 plot(recx{sp},recy{sp},'b')
 hold on
end
% scatter(x_corner_boundary(1,:),x_corner_boundary(2,:),'filled');

%% Particle variables
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
d_sp                    = zeros(spCount,2);                     % displacement
b_sp                    = zeros(spCount,2);                     % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
v_ssp                   = zeros(spCount,2);                     % velocty
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation
% Gradient deformation
for sp = 1:spCount
    F_sp{sp} = [1 0; 0 1];
end

% Data for analytical solutions
R_sp        = zeros(spCount,1);
Constant    = zeros(spCount,1);
alfap       = zeros(spCount,1);
alpha       = zeros(spCount,1);
theta       = zeros(spCount,1);     
beta        = zeros(spCount,1);    
b_r         = zeros(spCount,1);
b_angle     = zeros(spCount,1);
vel_ana     = zeros(spCount,2);
dis_ana     = zeros(spCount,2);
PK_ana      = zeros(spCount,3);
Cauchy_ana  = zeros(spCount,3);

xx=zeros(spCount,1);    vv=zeros(spCount,1);    ss=zeros(spCount,1);

%% Initial condition
for sp=1:spCount
R_sp(sp)                = sqrt((x_spo(sp,1)-center)^2 + (x_spo(sp,2)-center)^2);
theta(sp)               = atan2((x_spo(sp,2)-center),(x_spo(sp,1)-center));
Constant(sp)            = 1 - 32 * (R_sp(sp)-1)^2 + 256 * (R_sp(sp)-1)^4;
alpha(sp)               = A * sin(C*pi*t) * Constant(sp);
alfap(sp)               = A^2 * C * pi * Constant(sp);
v_ssp(sp,1)             = alfap(sp) * (-sin(alpha(sp))*(x_spo(sp,1)-center) - cos(alpha(sp))*(x_spo(sp,2)-center));
v_ssp(sp,2)             = alfap(sp) * (cos(alpha(sp))*(x_spo(sp,1)-center) - sin(alpha(sp))*(x_spo(sp,2)-center));
end

m_sp                    = psp * V_sp;                           % mass

%% Build velocity data for IMLS
velocity_data = cell(spCount,5);

for sp = 1:spCount  
    for i = 1:5
    R                       = sqrt((x_data_initial{sp,i}(1)-center)^2 + (x_data_initial{sp,i}(2)-center)^2);
    Constant1                = 1 - 32 * (R-1)^2 + 256 * (R-1)^4;
    alpha1                   = A * sin(C*pi*t) * Constant1;
    alfap1                   = A^2 * C * pi * Constant1;
    velocity_data{sp,i}(1,1)  = alfap1 * (-sin(alpha1)*(x_data_initial{sp,i}(1)-center) - cos(alpha1)*(x_data_initial{sp,i}(2)-center));
    velocity_data{sp,i}(2,1)  = alfap1 * (cos(alpha1)*(x_data_initial{sp,i}(1)-center) - sin(alpha1)*(x_data_initial{sp,i}(2)-center));
    end
end

%% Build velocity gradient
C_sp                    = cell(spCount,1);
 for sp = 1:spCount
        C_sp{sp} = zeros(2,2);
 end
%  
% start the algorithm
% video
timestep = 20;     % number of frame to save
r=timestep/5;      % number of frame per second video ~200s

writerObj2           = VideoWriter('Vortex2D.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:timestep
    ft              = ftime/timestep*tt;
%     ft=ftime;
 while t<ft+0.0001      
     Version
     resolution
     t
     
    % Methods of Manufactued solution
    % Compute the analytical solutions
    for sp = 1:spCount
        alpha(sp)   = A*sin(C*pi*(n-1)*dt)*Constant(sp);
        beta(sp)    = theta(sp) + alpha(sp);        
        b_r(sp)     = -R_sp(sp)*A^2*(cos(C*pi*(n-1)*dt))^2*C^2*pi^2*(4*R_sp(sp)-3)^4*(4*R_sp(sp)-5)^4+2/psp*2048*(4*R_sp(sp)-5)^2*(R_sp(sp)-1)^2*(4*R_sp(sp)-3)^2*A^2*(sin(C*pi*(n-1)*dt))^2*R_sp(sp)*Mu;
        b_angle(sp) = -A*sin(C*pi*(n-1)*dt)*(R_sp(sp)*C^2*pi^2*(4*R_sp(sp)-3)^2*(4*R_sp(sp)-5)^2+1/psp*64*Mu*(96*R_sp(sp)^3-240*R_sp(sp)^2+188*R_sp(sp)-45));
        b_sp(sp,1)  = b_r(sp) * cos(beta(sp)) - b_angle(sp) * sin(beta(sp));
        b_sp(sp,2)  = b_r(sp) * sin(beta(sp)) + b_angle(sp) * cos(beta(sp));     
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
    
        [gridmomentum,particlemomentum,k,color_node,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,...
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
    
    for sp = 1:spCount
      
        alpha(sp)   = A*sin(C*pi*n*dt)*Constant(sp);
        beta(sp)    = theta(sp) + alpha(sp); 
        
    % Analytical solution
    alfap(sp)       = A^2 * C * pi * Constant(sp);
    vel_ana(sp,1)   = alfap(sp) * cos(C*pi*n*dt) * (-sin(alpha(sp))*(x_spo(sp,1)-center) - cos(alpha(sp))*(x_spo(sp,2)-center));
    vel_ana(sp,2)   = alfap(sp) * cos(C*pi*n*dt) * (cos(alpha(sp))*(x_spo(sp,1)-center) - sin(alpha(sp))*(x_spo(sp,2)-center));
    
    dis_ana(sp,1)   = (A * (cos(alpha(sp))*(x_spo(sp,1)-center) - sin(alpha(sp))*(x_spo(sp,2)-center))) - (x_spo(sp,1)-center);
    dis_ana(sp,2)   = (A * (sin(alpha(sp))*(x_spo(sp,1)-center) + cos(alpha(sp))*(x_spo(sp,2)-center))) - (x_spo(sp,2)-center);
    
    PK_ana(sp,1) = 0;
    PK_ana(sp,2) = A^2 * R^2 * Mu * sin((C*pi*n*dt)^2) * (-64*(-1+R)+1024*(-1+R)^3)^2;
    PK_ana(sp,3) = A * R * Mu * sin(C*pi*n*dt) * (-64*(-1+R)+1024*(-1+R)^3);
    
    Cauchy_ana(sp,1) = sin(beta(sp)^2)*PK_ana(sp,2) - 2*sin(beta(sp))*cos(beta(sp))*PK_ana(sp,3);
    Cauchy_ana(sp,2) = cos(beta(sp)^2)*PK_ana(sp,2) + 2*sin(beta(sp))*cos(beta(sp))*PK_ana(sp,3);
    Cauchy_ana(sp,3) = -sin(beta(sp))*cos(beta(sp))*PK_ana(sp,2) + ((cos(beta(sp)))^2-(sin(beta(sp)))^2) * PK_ana(sp,3);
    end
    
    % Update time and step 
     t = t+dt;
     n = n+1;

 end

%     deviation_V(sp,:) = vel_ana(sp,:)-v_ssp(sp,:);
%     v=v + norm(deviation_V(sp,:))^2;
%     vv = norm(deviation_V(sp,:));
%     
%     deviation_S(sp,:) = Cauchy_ana(sp,:)-s_sp(sp,:);
%     s=s + norm(deviation_S(sp,:))^2;
%     ss = norm(deviation_S(sp,:));
%     end
%     
%     perror_X = sqrt(x/spCount); perror_V = sqrt(v/spCount); perror_S = sqrt(s/spCount);   
%     results = [perror_X perror_V perror_S];
% 
%     
%    %% Plot the result
% %     StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,color_node,v_ssp,spCount,spElems_boundary,cellCount,nodeCount,r1_sp,r2_sp,p_sc,LOCC);  
%     %% Plot particle
% 
% %     velocity = zeros(spCount,1);
% %     for sp=1:spCount
% %     displacement(sp) = sqrt(d_sp(sp,1)^2+d_sp(sp,2)^2);
% %     velocity(sp) = sqrt(v_ssp(sp,1)^2+v_ssp(sp,2)^2);
% %     end
% x_corner1 = x_sp - r1_sp - r2_sp;
% x_corner2 = x_sp + r1_sp - r2_sp;
% x_corner3 = x_sp + r1_sp + r2_sp;
% x_corner4 = x_sp - r1_sp + r2_sp;
% 
%     StressProfile1=figure;
%     sz = 100*le(1);
%     color = xx;
%     scatter(x_sp(:,1),x_sp(:,2),sz,color,'filled');
%     hold on
%     %% Plot particle domain
%     for sp=1:spCount
%     % %     Qxx{sp} = [x_sp(sp,1) x_sp(sp,1)+r1_sp(sp,1)];
%     % %     Qxy{sp} = [x_sp(sp,2) x_sp(sp,2)+r1_sp(sp,2)];
%     % %     Qyx{sp} = [x_sp(sp,1) x_sp(sp,1)+r2_sp(sp,1)];
%     % %     Qyy{sp} = [x_sp(sp,2) x_sp(sp,2)+r2_sp(sp,2)];
%     % %     plot (Qxx{sp},Qxy{sp},'r')
%     % %     hold on 
%     % %     plot (Qyx{sp},Qyy{sp},'r')
%     % %     hold on
%     %     
%         x_cor = [x_corner1(sp,1) x_corner2(sp,1) x_corner3(sp,1) x_corner4(sp,1) x_corner1(sp,1)];
%         y_cor = [x_corner1(sp,2) x_corner2(sp,2) x_corner3(sp,2) x_corner4(sp,2) x_corner1(sp,2)];
%         plot(x_cor,y_cor,'b')
%         hold on
%     end
%     grid on
% 
%     %% Fornat
%     axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
%     % axis([1.25,3.75,1.25,3.75]);
%     % axis([0.9,1.5,1.5,3]);
%     set(gca,'xtick',[0:le(1):max(LOC(:,1))]);
%     set(gca,'ytick',[0:le(2):max(LOC(:,2))]);
%     h=colorbar;
%     colormap(jet(256))
%     % set(h, 'ytick', [0:0.2:1.2]);
% %         caxis([0 0.1]);
% 
% frame2 = getframe(StressProfile1);
%     writeVideo(writerObj2,frame2);
% %     close(StressProfile1);
    end
%     close(writerObj2);

     %% Compute error
    x=0;    v=0;    s=0;
    deviation_X = zeros(spCount,2);    deviation_V = zeros(spCount,2);    deviation_S = zeros(spCount,3);
    
    % Loop all particles to store the RMS
    for sp=1:spCount
    deviation_X(sp,:) = dis_ana(sp,:)-d_sp(sp,:);
    x=x + norm(deviation_X(sp,:))^2;
    
    deviation_V(sp,:) = vel_ana(sp,:)-v_ssp(sp,:);
    v=v + norm(deviation_V(sp,:))^2;
    
    deviation_S(sp,:) = Cauchy_ana(sp,:)-s_sp(sp,:);
    s=s + norm(deviation_S(sp,:))^2;
    end
    
    perror_X = sqrt(x/spCount); perror_V = sqrt(v/spCount); perror_S = sqrt(s/spCount);   
    results = [perror_X perror_V perror_S];
