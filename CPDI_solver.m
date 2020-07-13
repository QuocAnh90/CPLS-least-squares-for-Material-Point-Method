function[gridmomentum,particlemomentum,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp,nvelo_si] = CPDI_solver(CModel,CModel_parameter,...
    nodeCount,spCount,cellCount...
    ,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
    nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp,gridmomentum,particlemomentum)
     
%% Store particles into cell
[N,dN,CONNECT,spElems,mspoints,NODES] = Compute_Interpolator_CPDI(spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp,V_sp);
% N: Shape function
% dN: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

 %% Mapping from particle to nodes
% [nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si]=Interpolate_Particle_To_Grid(NODES,nodeCount,CONNECT,le,N,dN,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp);
% nmass_si: nodal mass
% nmomentum_si: nodal momentum
% niforce_si: nodal internal force
% neforce_si: nodal external force
% traction_si: nodal traction

%% Interpolation from particle to grid task
% Node variables
nmass_si                = zeros(nodeCount,1);                   % Nodal Mass
nmomentum_si            = zeros(nodeCount,2);                   % Nodal Momentum
niforce_si              = zeros(nodeCount,2);                   % Nodal Internal force
neforce_si              = zeros(nodeCount,2);                   % Nodal External force
traction_si             = zeros(nodeCount,2);                   % Nodal Traction
particle = zeros(spCount,2);
% node = zeros(spCount,2);
% test = zeros(spCount,2);

 for p=1:spCount
 % Build stress tensor
 SSP = [s_sp(p,1) s_sp(p,3);s_sp(p,3) s_sp(p,2)];
 
 for j=1:NODES(p)
     npid                           = CONNECT{p}(j);
     
          if N{p}(j)==0
         continue
          end
     
 % Mass
 nmass_si(npid)            = nmass_si(npid) + m_sp(p)*N{p}(j);
 
 % Momentum
 nmomentum_si(npid,:)      = nmomentum_si(npid,:) + m_sp(p)*v_ssp(p,:)*N{p}(j);
%  particle(p,:)  = particle(p,:) + m_sp(p)*v_ssp(p,:)*N{p}(j); 
 
 % Internal forces
niforce_si(npid,:)         = niforce_si(npid,:) - (V_sp(p)*SSP*dN{p}(:,j))';

 % External forces
neforce_si(npid,:)         = neforce_si(npid,:) + b_sp(p,:)*m_sp(p)*N{p}(j);

% Traction
%  traction_si(npid,:)       = traction_si(npid,:) + V_sp(sp)*ptraction_sp(sp,:)*N{sp}(j)/le(1,1)/le(1,2);
 traction_si(npid,1)       = traction_si(npid,1) + V_sp(p)*ptraction_sp(p,1)*N{p}(j)/le(1);
 traction_si(npid,2)       = traction_si(npid,2) + V_sp(p)*ptraction_sp(p,2)*N{p}(j)/le(2); 
 end 
%   node(p,:)      = m_sp(p)*v_ssp(p,:);
% test(p,:) = sum(nmomentum_si(:,1)) - sum(node(:,1));
 end
 
  
momentumgrid = 0;
% for n = 1:nodeCount
%     momentumgrid = momentumgrid + abs((nmomentum_si(n,1)));
% end
momentumparticle = 0;
% for p = 1:spCount
%     momentumparticle = momentumparticle + norm(abs(node(p,:)));
%     momentumgrid = momentumgrid + norm(abs((particle(p,:))));
% end
% 
% gridmomentum = [gridmomentum;momentumgrid];
% particlemomentum = [particlemomentum;momentumparticle];


%% Update momentum
% Update force and momentum
 nforce_si      = niforce_si + neforce_si + traction_si;
%  nmomentum_old = nmomentum_si;
 nmomentum_si   = nmomentum_si + nforce_si*dt;
 
% Boundary condition
[nforce_si]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nforce_si); % Boundary condition for nodal force
[nmomentum_si]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nmomentum_si); % Boundary condition for nodal force

%     nforce_si = (nmomentum_si - nmomentum_old)/dt;

%% Update solid particle velocity and position
[v_ssp,x_sp,d_sp] = Update_Particle_Position(NODES,dt,CONNECT,N,spCount,nmass_si,nforce_si,nmomentum_si,x_spo,v_ssp,x_sp);
% velocity particle: v_ssp
% position particle: x_sp
% displacement particle: d_sp

%% Mapping nodal velocity back to node
% [nvelo_si] = Interpolate_velocity_back(NODES,nodeCount,spCount,CONNECT,m_sp,v_ssp,N,nmass_si);
nvelo_si = zeros(nodeCount,2);
for n = 1:nodeCount
    if nmass_si(n) >1e-12
    nvelo_si(n,:) = nmomentum_si(n,:)/nmass_si(n);
    end
end

% Boundary condition
[nvelo_si] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si); % Boundary condition for nodal force
% [~,dN,CONNECT,~,mspoints,NODES] = Compute_Interpolator_CPDI(spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp,V_sp);

%% Update effective stress
[F_sp,V_sp,s_sp,p_sp,L_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES,dt,cellCount,mspoints,CONNECT,nvelo_si,dN,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);

%% Update the topology of particles
[r1_sp,r2_sp] = Update_topology(spCount,F_sp,r1_sp,r10_sp,r2_sp,r20_sp);