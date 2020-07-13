% function results = Vibration2D(Version,resolution)

clear all
close all
addpath('SubFunctions');
% Unit
% Newton - seconds - metre

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
Version = 'CPDI';

%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
CModel = 'Neo_Hookean_Elastic';

% Amplitude of displacement
A = 0.1;

%% Material porperties
E                       = 10000000          ;                % Young modulus of solid
psp                     = 1000.0            ;               % solid density
nu                      = 0.3           ;                   % Poison ratio
g                       = 0.0            ;                  % gravity acceleration
Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);
CModel_parameter = [E,nu];

%% Structured Grid input
resolution              = 8;
NN(1)                   = resolution + 5;                  % number of nodes in X direction
NN(2)                   = resolution + 5;                  % number of nodes in Y direction
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
ftime                   = 0.02;
dt                      = 0.4*le(1)/C;
ndt                     = round(ftime/dt) +1;
t                       = 0;

% Plot
spElems_boundary = 0;
p_sc = zeros(cellCount,1);
color_node = zeros(nodeCount,1);
gridmomentum=[];
particlemomentum=[];

%% Boundary nodes
% Boundary coordination
x_min = 2*le(1);
x_max = 1+2*le(1);
y_min = 2*le(2);
y_max = 1+2*le(2);
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

%% Particle generation
particle_per_cell       = 4;
spCount                 = resolution*resolution*particle_per_cell;
lp(1)                   = le(1)/sqrt(particle_per_cell);                                 % size of particle in X direction
lp(2)                   = le(2)/sqrt(particle_per_cell);                                % size of particle in Y direction
x_sp                    = zeros(spCount,2);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:resolution*sqrt(particle_per_cell)
        for j=1:resolution*sqrt(particle_per_cell)
            x_sp(sp,1:2)= [2*le(2)+0.5*lp(1,1)+(j-1)*lp(1,1) 2*le(2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
            sp=sp+1;
        end
    end  
end
% x_sp: Vector, position of MPs
% spCount: total number of MPs

%% Plot initial condition
% initial_figure = Plot_Initial(x_sp,LOC,le);

%% Particle variables
dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
d_sp                    = zeros(spCount,2);                     % displacement
b_sp                    = zeros(spCount,2);                     % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
v_ssp1                  = zeros(spCount,2);                     % velocty
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
boundary_corner = [];

for sp = 1:spCount
    x_data{sp,1} = x_sp(sp,:)';
    x_data{sp,2} = x_corner1(sp,:)';
    x_data{sp,3} = x_corner2(sp,:)';
    x_data{sp,4} = x_corner3(sp,:)';
    x_data{sp,5} = x_corner4(sp,:)';
    
   % Determine the object boundary corners
    for i = 1:5
        X1 = x_data{sp,i}(1)-2*le(1); Y1 = x_data{sp,i}(2)-2*le(1);   
        
        if X1 < 0.1*le(1)
            boundary_corner = [boundary_corner [sp;i]];
        end
        
        if X1 > 1-0.1*le(1)
            boundary_corner = [boundary_corner [sp;i]];
        end 
        
        if Y1 < 0.1*le(2)
            boundary_corner = [boundary_corner [sp;i]];
        end
        
        if Y1 > 1-0.1*le(2)
            boundary_corner = [boundary_corner [sp;i]];
        end 
        
        
    end
end

x_data_initial = x_data;

% Take the unique value of boundary nodes
boundary_corner = boundary_corner';
boundary_corner = unique(boundary_corner,'rows');
boundary_corner = boundary_corner';

%% Generate the vibration
for sp=1:spCount
    v_ssp(sp,1) = A*C*pi*sin(pi*(x_spo(sp,1)-2*le(1)))*cos(C*pi*t);
    v_ssp(sp,2) = A*C*pi*sin(pi*(x_spo(sp,2)-2*le(2)))*cos(C*pi*t+pi);
end

%% Build velocity data for IMLS
velocity_data = cell(spCount,5);
ux_data = cell(spCount,5);
uy_data = cell(spCount,5);
Fxx_data = cell(spCount,5);
Fyy_data = cell(spCount,5);

for sp = 1:spCount  
    velocity_data{sp,1} = [A*C*pi*sin(pi*(x_spo(sp,1)-2*le(1)));A*C*pi*sin(pi*(x_spo(sp,2)-2*le(1)))];
    velocity_data{sp,2} = [A*C*pi*sin(pi*(x_corner1(sp,1)-2*le(1)));A*C*pi*sin(pi*(x_corner1(sp,2)-2*le(1)))];
    velocity_data{sp,3} = [A*C*pi*sin(pi*(x_corner2(sp,1)-2*le(1)));A*C*pi*sin(pi*(x_corner2(sp,2)-2*le(1)))];
    velocity_data{sp,4} = [A*C*pi*sin(pi*(x_corner3(sp,1)-2*le(1)));A*C*pi*sin(pi*(x_corner3(sp,2)-2*le(1)))];
    velocity_data{sp,5} = [A*C*pi*sin(pi*(x_corner4(sp,1)-2*le(1)));A*C*pi*sin(pi*(x_corner4(sp,2)-2*le(1)))];  
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
 while t<ft+0.0001      
     Version
     resolution
     t

      % Methods of Manufactued solution
     % Compute the analytical solutions
     for sp=1:spCount      
        % Displacement
        ux(sp) = A*sin(pi*(x_spo(sp,1)-2*le(1)))*sin(C*pi*t);
        uy(sp) = A*sin(pi*(x_spo(sp,2)-2*le(2)))*sin(C*pi*t+pi);   
        % Deformation gradient
        Fxx(sp) = 1 + A*pi*cos(pi*(x_spo(sp,1)-2*le(1)))*sin(C*pi*t);
        Fyy(sp) = 1 + A*pi*cos(pi*(x_spo(sp,2)-2*le(2)))*sin(C*pi*t+pi);
        % Body force applied to the particles
        b_sp(sp,1) = pi^2*ux(sp)/psp*(Lambda/Fxx(sp)^2*(1-log(Fxx(sp)*Fyy(sp)))+Mu/Fxx(sp)^2*(Fxx(sp)^2+1)-E);
        b_sp(sp,2) = pi^2*uy(sp)/psp*(Lambda/Fyy(sp)^2*(1-log(Fxx(sp)*Fyy(sp)))+Mu/Fyy(sp)^2*(Fyy(sp)^2+1)-E); 
        
        % Build body force for IMLS_CPDI
        for i = 1:5
        ux_data{sp,i} = A*sin(pi*(x_data_initial{sp,i}(1)-2*le(1)))*sin(C*pi*t);
        uy_data{sp,i} = A*sin(pi*(x_data_initial{sp,i}(2)-2*le(2)))*sin(C*pi*t+pi); 
        Fxx_data{sp,i} = 1 + A*pi*cos(pi*(x_data_initial{sp,i}(1)-2*le(1)))*sin(C*pi*t);  
        Fyy_data{sp,i} = 1 + A*pi*cos(pi*(x_data_initial{sp,i}(2)-2*le(2)))*sin(C*pi*t+pi);
        end
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
        [gridmomentum,particlemomentum,color_node,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,...
        r1_sp,r2_sp,L_sp,velocity_data,x_data]...
        = IMLS_CPDI_solver_test(gridmomentum,particlemomentum,...
        velocity_data,x_data,...
        CModel,CModel_parameter,...
        nodeCount,spCount,cellCount...
        ,x_sp,x_spo,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
    
     end
     
 % Update time and step 
 t = t+dt;
 
     % Methods of Manufactued solution
     % Compute the analytical solutions
     for sp=1:spCount      
        % Displacement
        ux(sp) = A*sin(pi*(x_spo(sp,1)-2*le(1)))*sin(C*pi*t);
        uy(sp) = A*sin(pi*(x_spo(sp,2)-2*le(2)))*sin(C*pi*t+pi);   
        % Velocity
        vx(sp) = A*C*pi*sin(pi*(x_spo(sp,1)-2*le(1)))*cos(C*pi*t);
        vy(sp) = A*C*pi*sin(pi*(x_spo(sp,2)-2*le(2)))*cos(C*pi*t+pi); 
        % Deformation gradient
        Fxx(sp) = 1 + A*pi*cos(pi*(x_spo(sp,1)-2*le(1)))*sin(C*pi*t);
        Fyy(sp) = 1 + A*pi*cos(pi*(x_spo(sp,2)-2*le(2)))*sin(C*pi*t+pi);
        % Stress
        Fa{sp} = [Fxx(sp) 0;0 Fyy(sp)];
        Ja(sp) = det(Fa{sp});
         
        [S{sp}]=Neo_Hookean_elastic(CModel_parameter, Fa{sp},Ja(sp));
        sx(sp) = S{sp}(1);
        sy(sp) = S{sp}(2);

        % Store analytical solutions 
        dis_ana(sp,:) =  [ux(sp) uy(sp)];
        vel_ana(sp,:) =  [vx(sp) vy(sp)];
        stress_ana(sp,:) =  [sx(sp) sy(sp) 0];
     end
 end

 %% Plot the result
StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,color_node,v_ssp,spCount,spElems_boundary,cellCount,nodeCount,r1_sp,r2_sp,p_sc,LOCC);  
    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    close(writerObj2);
%     
     %% Compute error
    x=0;    v=0;    s=0;
    deviation_X = zeros(spCount,2);    deviation_V = zeros(spCount,2);    deviation_S = zeros(spCount,3);
    
    % Loop all particles to store the RMS
    for sp=1:spCount
    deviation_X(sp,:) = dis_ana(sp,:)-d_sp(sp,:);
    x=x + norm(deviation_X(sp,:))^2;
    
    deviation_V(sp,:) = vel_ana(sp,:)-v_ssp(sp,:);
    v=v + norm(deviation_V(sp,:))^2;
    
    deviation_S(sp,:) = stress_ana(sp,:)-s_sp(sp,:);
    s=s + norm(deviation_S(sp,:))^2;
    end
    
    perror_X = sqrt(x/spCount); perror_V = sqrt(v/spCount); perror_S = sqrt(s/spCount);   
    results = [perror_X perror_V perror_S];