function results = Vibration2D(Version,resolution)
% clear all
% close all

% Unit
% Newton - seconds - metre
addpath('SubFunctions');

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'MPM'
% B spline MPM: 'BSMPM'
% Moving least square MPM: 'MLS_BSMPM'
% Version = 'MLS_BSMPM';

%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
CModel = 'Neo_Hookean_Elastic';

% Amplitude of displacement
A = 0.1;
particle_index = 2;

%% Material porperties
E                       = 10000000          ;                % Young modulus of solid
psp                     = 1000.0            ;               % solid density
nu                      = 0.3           ;                   % Poison ratio
g                       = 0.0            ;                  % gravity acceleration
Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);
CModel_parameter = [E,nu];

%% Structured Grid input
% resolution              = 8;
NN(1)                   = resolution + 5;                               % number of nodes in X direction
NN(2)                   = resolution + 5;                               % number of nodes in Y direction
le(1)                   = 1/resolution;                              % size of element in X direction
le(2)                   = 1/resolution;                            % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

% Generate the Gauss point
GaussCount = 4;
x_gauss = zeros(cellCount, GaussCount);
y_gauss = zeros(cellCount, GaussCount);

xc_gauss = zeros(cellCount, 9*GaussCount);
yc_gauss = zeros(cellCount, 9*GaussCount);
omega_gauss = zeros(cellCount, GaussCount);

    for c=1:cellCount
    
    if GaussCount==1
        n_g = 1;
    else
        n_g = GaussCount/2;
    end
    
    x0 = LOCC(c,1)-0.5*le(1);               
    x1 = LOCC(c,1)+0.5*le(1);    
    [x_gauss_temp, x_omega_gauss] = lgwt(n_g,x0,x1);
    
    if GaussCount==4
    x_gauss(c,:) = [x_gauss_temp(1), x_gauss_temp(2), x_gauss_temp(2), x_gauss_temp(1)];   
    elseif GaussCount==1
    x_gauss(c,:) = [x_gauss_temp(1)];
    end
    
    y0 = LOCC(c,2)-0.5*le(2);               
    y1 = LOCC(c,2)+0.5*le(2);   
    [y_gauss_temp, y_omega_gauss] = lgwt(n_g,y0,y1);
    
    if GaussCount==4
    y_gauss(c,:) = [y_gauss_temp(1), y_gauss_temp(1), y_gauss_temp(2), y_gauss_temp(2)];
    elseif GaussCount==1
    y_gauss(c,:) = [y_gauss_temp(1)];
    end
    
    omega_gauss_temp = x_omega_gauss*y_omega_gauss';
    omega_gauss(c,:) = reshape(omega_gauss_temp, 1, GaussCount);
    end
    

%% Time
C                       = sqrt(E/psp);
ftime                   = 0.02;
dt                      = 0.00001;
ndt                     = round(ftime/dt) +1;
t                       = 0;

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
particle_per_cell       = 3;
spCount                 = 3*resolution*resolution*particle_per_cell;
lp(1)                   = le(1)/(particle_per_cell);              % size of particle in X direction
lp(2)                   = le(2)/(particle_per_cell);              % size of particle in Y direction
x_sp                    = zeros(spCount,2);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:resolution*particle_per_cell
        for j=1:resolution*particle_per_cell
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

% Generate the vibration
for sp=1:spCount
    v_ssp(sp,1) = A*C*pi*sin(pi*(x_spo(sp,1)-2*le(1)))*cos(C*pi*t);
    v_ssp(sp,2) = A*C*pi*sin(pi*(x_spo(sp,2)-2*le(2)))*cos(C*pi*t+pi);
end

% Data
sx1=[];sy1=[];sx2=[];sy2=[];
ux1=[];uy1=[];
time=[];
error_X=[];error_V=[];error_S=[];
perror_X=0;perror_V=0;perror_S=0;
dx1=[];dy1=[];
Fxx1=[];Fxx2=[];Fyy1=[];Fyy2=[];
bxx=[];byy=[];

% % start the algorithm
% % video
% timestep = 100;     % number of frame to save
% r=timestep/10;      % number of frame per second video ~200s
% 
% writerObj2           = VideoWriter('aligned1Dx.avi');
% writerObj2.FrameRate = r;    % number of frame per second
% open(writerObj2);
% 
%     for tt = 1:timestep
%     ft              = ftime/timestep*tt;
    ft=ftime;
 while t<ft+0.0000001      
     t

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
        
        % Store analytical solutions5
        dis_ana(sp,:) =  [ux(sp) uy(sp)];
        vel_ana(sp,:) =  [vx(sp) vy(sp)];
        stress_ana(sp,:) =  [sx(sp) sy(sp) 0];
        
       % Body force applied to the particles
        b_sp(sp,1) = pi^2*ux(sp)/psp*(Lambda/Fxx(sp)^2*(1-log(Fxx(sp)*Fyy(sp)))+Mu/Fxx(sp)^2*(Fxx(sp)^2+1)-E);
        b_sp(sp,2) = pi^2*uy(sp)/psp*(Lambda/Fyy(sp)^2*(1-log(Fxx(sp)*Fyy(sp)))+Mu/Fyy(sp)^2*(Fyy(sp)^2+1)-E);
     end
     
     switch Version
         case 'MPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MPM_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,lp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);
        
        case 'BSMPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = BSMPM_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);
        
        case 'MLS_BSMPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MLS_BSMPM_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOCC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        F_sp,V_spo,dt,GaussCount,x_gauss,y_gauss,omega_gauss);

     end
       
    time = [time t];
    % Store for 1st particle     
    % Analytical solutions for particle 1
    ux1 = [ux1 ux(particle_index)];    sx1 = [sx1 sx(particle_index)];
    Fxx1 = [Fxx1 Fxx(particle_index)];    bxx = [bxx b_sp(particle_index,1)];
    
    uy1 = [uy1 uy(particle_index)];    sy1 = [sy1 sy(particle_index)];
    Fyy1 = [Fyy1 Fyy(particle_index)];    byy = [byy b_sp(particle_index,2)];
    
    % Numerical solutions for particle 1
    sx2 = [sx2 s_sp(particle_index,1)];    dx1 = [dx1 d_sp(particle_index,1)];
    Fxx2 = [Fxx2 F_sp{particle_index}(1,1)];
    
    sy2 = [sy2 s_sp(particle_index,2)];    dy1 = [dy1 d_sp(particle_index,2)];
    Fyy2 = [Fyy2 F_sp{particle_index}(2,2)];
    
 % Update time and step 
 t = t+dt;
 end

%  %% Plot the result
% StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,0,v_ssp,spCount,r1_sp,r2_sp);
%     frame2 = getframe(StressProfile1);
%     writeVideo(writerObj2,frame2);
%     end
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
    
    deviation_S(sp,:) = stress_ana(sp,:)-s_sp(sp,:);
    s=s + norm(deviation_S(sp,:))^2;
    end
    
    perror_X = sqrt(x/spCount); perror_V = sqrt(v/spCount); perror_S = sqrt(s/spCount);   
    results = [perror_X perror_V perror_S];
    
%     figure
%     plot(time, ux1, time, dx1);
%     title('ux');
%     legend('analytical','numerical')
% 
%     figure 
%     plot(time,sx1,time,sx2);
%     title('sx');
%     legend('analytical','numerical')
%     
%     figure 
%     plot(time,Fxx1,time,Fxx2);
%     title('Fx');
%     legend('analytical','numerical')
%     
%     figure
%     plot(time,bxx);
%     title('bx');
%     legend('analytical','numerical')
%     
%       figure
%     plot(time, uy1, time, dy1);
%     title('uy');
%     legend('analytical','numerical')
% 
%     figure 
%     plot(time,sy1,time,sy2);
%     title('sy');
%     legend('analytical','numerical')
%     
%     figure 
%     plot(time,Fyy1,time,Fyy2);
%     title('Fy');
%     legend('analytical','numerical')
%     
%     figure
%     plot(time,byy);
%     title('by');
%     legend('analytical','numerical')
%     
% figure
% plot(x_sp(1:resolution*particle_per_cell,1),stress_ana(1:resolution*particle_per_cell,1),x_sp(1:resolution*particle_per_cell,1),s_sp(1:resolution*particle_per_cell,1))
% legend('analytical','numerical')
% xlabel('position')
% ylabel('stress')
% 
% figure
% plot(x_sp(1:resolution*particle_per_cell,1),vel_ana(1:resolution*particle_per_cell,1),x_sp(1:resolution*particle_per_cell,1),v_ssp(1:resolution*particle_per_cell,1))
% legend('analytical','numerical')
% xlabel('position')
% ylabel('velocity')
% 
% figure
% plot(x_sp(1:resolution*particle_per_cell,1),dis_ana(1:resolution*particle_per_cell,1),x_sp(1:resolution*particle_per_cell,1),d_sp(1:resolution*particle_per_cell,1))
% legend('analytical','numerical')
% xlabel('position')
% ylabel('displacement')
