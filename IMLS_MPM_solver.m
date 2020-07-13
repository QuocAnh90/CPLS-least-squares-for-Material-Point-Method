function[v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = IMLS_MPM_solver(CModel,CModel_parameter,...
    nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,LOCC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
    nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt)

%% Store particles into cell
[N,dN,N_MLS,N_g,N_c,dN_c,CONNECT,CONNECT_G,CONNECT_C,CONNECT_N...
,spElems,mspoints,NODES,GAUSS,NODES_G,NODES_N] = Compute_Interpolator_IMLS_MPM(spCount,cellCount,nodeCount,x_sp,le,NN,LOC,LOCC);
% N: Shape function
% dN: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

%% Interpolation from particle to grid task
[nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si]=Interpolate_Particle_To_Grid_MLS(NODES,NODES_N,GAUSS,NODES_G,cellCount,nodeCount,...
    CONNECT,CONNECT_G,CONNECT_C,CONNECT_N,...
    le,N,dN,N_g,N_c,dN_c,N_MLS,...
    spCount,m_sp,b_sp,V_sp,p_sp,ptraction_sp,v_ssp,s_sp,mspoints);

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
% [~,dN,N_MLS,~,~,~,CONNECT,~,~,CONNECT_N...
% ,~,mspoints,NODES,~,~,NODES_N] = Compute_Interpolator_IMLS_MPM(spCount,cellCount,nodeCount,x_sp,le,NN,LOC,LOCC);

% [nvelo_si] = Interpolate_velocity_back_MLS(NODES,NODES_N,nodeCount,spCount,CONNECT,CONNECT_N,v_ssp,N_MLS);
nvelo_si = zeros(nodeCount,2);
for n = 1:nodeCount
    if nmass_si(n) == 0
        continue
    end
    nvelo_si(n,:) = nmomentum_si(n,:)/nmass_si(n);
end

% Boundary condition
[nvelo_si] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si); % Boundary condition for nodal force

%% Update effective stress
[F_sp,V_sp,s_sp,p_sp,L_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES,dt,cellCount,mspoints,CONNECT,nvelo_si,dN,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);
