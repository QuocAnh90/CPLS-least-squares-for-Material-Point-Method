function[gridmomentum,particlemomentum,k,color_node,v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,...
        r1_sp,r2_sp,L_sp,velocity_data,x_data,nvelo_si]...
        = IMLS_CPDI_solver_test(gridmomentum,particlemomentum,...
        velocity_data,x_data,...
        CModel,CModel_parameter,...
        nodeCount,spCount,cellCount...
        ,x_sp,x_spo,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp)
    
%% Store particles into cell
    [k,node_boundary,color_cell,color_node,N_MLS,...
    CONNECT_MLS]...
    = Compute_Interpolator_IMLS_CPDIv2(x_data,...
    nodeCount,spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp);

[N,dN,CONNECT,spElems,mspoints,NODES] = Compute_Interpolator_CPDI...
    (spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp,V_sp);

%% Mapping from particle to cell
% Node variables
nmass_si                = zeros(nodeCount,1);                   % Nodal Mass
niforce_si              = zeros(nodeCount,2);                   % Nodal Internal force
neforce_si              = zeros(nodeCount,2);                   % Nodal External force
traction_si             = zeros(nodeCount,2);                   % Nodal Traction
nmomentum_si            = zeros(nodeCount,2);                   % Nodal Momentum

node = zeros(spCount,2);
particle = zeros(spCount,2);

 for p=1:spCount
 % Build stress tensor
 SSP = [s_sp(p,1) s_sp(p,3);s_sp(p,3) s_sp(p,2)];
 
 for j=1:NODES(p)
     npid                           = CONNECT{p}(j);
     
          if N{p}(j)==1
         continue
          end
     
 % Mass
 nmass_si(npid)            = nmass_si(npid) + m_sp(p)*N{p}(j);
 
 % Internal forces
niforce_si(npid,:)         = niforce_si(npid,:) - (V_sp(p)*SSP*dN{p}(:,j))';

 % External forces
neforce_si(npid,:)         = neforce_si(npid,:) + b_sp(p,:)*m_sp(p)*N{p}(j);

 % Traction
%  traction_si(npid,:)       = traction_si(npid,:) + V_sp(sp)*ptraction_sp(sp,:)*N{sp}(j)/le(1,1)/le(1,2);
 traction_si(npid,1)       = traction_si(npid,1) + V_sp(p)*ptraction_sp(p,1)*N{p}(j)/le(1);
 traction_si(npid,2)       = traction_si(npid,2) + V_sp(p)*ptraction_sp(p,2)*N{p}(j)/le(2); 
 end 
 
%  node(p,:)      = m_sp(p)*v_ssp(p,:);
 end
 

%% Mapping velocity from particle to nodes
nvelo_si                = zeros(nodeCount,2);                   % Nodal Velocity

 for p = 1 : spCount        % Loop particle
     for c = 1:5            % Loop data including particle and corners
         for j= 1:16
             npid = CONNECT_MLS{p,c}(j);
             
             % Nodal Velocity
              nvelo_si(npid,:) = nvelo_si(npid,:) + velocity_data{p,c}' * N_MLS{p,c}(j);
%                particle(p,:)  = particle(p,:) + velocity_data{p,c}' * N_MLS{p,c}(j);

         end
     end
 end
 
 % Compute momentum
 for n = 1:nodeCount
%       if ismember(n,node_boundary)==0
      nmomentum_si(n,:)         = nvelo_si(n,:) * nmass_si(n);
%       end
 end
 
  momentumgrid = 0;
% for n = 1:nodeCount
%     momentumgrid = momentumgrid + norm(abs(nmomentum_si(n,:)));
% end
momentumparticle = 0;
% for p = 1:spCount
%     momentumparticle = momentumparticle + norm(abs(v_ssp(p,:).*m_sp(p)));
% end

% gridmomentum = [gridmomentum;momentumgrid];
% particlemomentum = [particlemomentum;momentumparticle];

%% Update momentum
% Update force and momentum
 nforce_si      = niforce_si + neforce_si + traction_si;
 nmomentum_si   = nmomentum_si + nforce_si*dt;
 
% Boundary condition
[nforce_si]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nforce_si); % Boundary condition for nodal force
[nmomentum_si]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nmomentum_si); % Boundary condition for nodal force

%% Update solid particle velocity and position
[v_ssp,x_sp,d_sp] = Update_Particle_Position(NODES,dt,CONNECT,N,spCount,nmass_si,nforce_si,nmomentum_si,x_spo,v_ssp,x_sp);
% velocity particle: v_ssp
% position particle: x_sp
% displacement particle: d_sp

%% Mapping nodal velocity back to node
nvelo_si = zeros(nodeCount,2);
for n = 1:nodeCount
    if nmass_si(n) >1e-12
    nvelo_si(n,:) = nmomentum_si(n,:)/nmass_si(n);
    end
end
 
% Boundary condition
[nvelo_si] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si); % Boundary condition for nodal force

%% Update effective stress
[F_sp,V_sp,s_sp,p_sp,L_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES,dt,cellCount,mspoints,CONNECT,nvelo_si,dN,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);

%% Update the topology of particles
[r1_sp,r2_sp] = Update_topology(spCount,F_sp,r1_sp,r10_sp,r2_sp,r20_sp);

%% Update velocity field at the corner
[velocity_data,x_data] = Gradient(spCount,x_sp,r1_sp,r2_sp,v_ssp,L_sp);
 