function[v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MLS_BSMPM_solvertest1(CModel,CModel_parameter,...
    nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOCC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
    F_sp,V_spo,dt,GaussCount,x_gauss,y_gauss,omega_gauss)

%% Store particles into cell
choice =2;
r = 3;
GaussCount = 2;
lump = 0;

if choice ==1
%% Quadratic Bspline
     deg = 2;
     Nelemsx = (NN(1)) - 4 + 1;
     XiCount = 1/le(1)+5; 
     Xi = zeros(XiCount,1);
     Xi(1:deg+1,1) = 2*le(1);
     for i=1:1/le(1)
     Xi(deg+1+i,1) =  Xi(deg+i,1)+le(1);
     end
     Xi(XiCount-1:XiCount) = Xi(XiCount-2);
     
    spElems         = zeros(spCount,1);
    CONNECT_TEMPX   = zeros(spCount,deg+1);
     for p = 1:spCount
     % compute vector store index elements 
     spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));

     % Knot vector index
     CONNECT_TEMPX(p,1) = floor(x_sp(p,1)/le(1)+1) - 2;
     CONNECT_TEMPX(p,2) = floor(x_sp(p,1)/le(1)+1) - 1;
     CONNECT_TEMPX(p,3) = floor(x_sp(p,1)/le(1)+1) - 0;
     CONNECTX{p}        = [CONNECT_TEMPX(p,1) CONNECT_TEMPX(p,2) CONNECT_TEMPX(p,3)];
     [N_localx , dN_localx] = test_new_bspline(Xi,deg,x_sp(p,1));
     
     for i = 1:deg+1
         Nx{p}(i) = N_localx(CONNECT_TEMPX(p,i));
         dNx{p}(1,i) = dN_localx(CONNECT_TEMPX(p,i));
     end  
     end

      %% Compute bspline shape function from centroids to knot
    active_elements = unique(spElems);
    ElemsCount = length(active_elements);
    CONNECT_TEMPXC   = zeros(ElemsCount,3);

    for c = 1:ElemsCount
    cpid = active_elements(c);
    
    % Knot vector index for centroids 
     CONNECT_TEMPXC(c,1) = floor(LOCC(cpid,1)/le(1)+1) - 2;
     CONNECT_TEMPXC(c,2) = floor(LOCC(cpid,1)/le(1)+1) - 1;
     CONNECT_TEMPXC(c,3) = floor(LOCC(cpid,1)/le(1)+1) + 0;
     CONNECTXC{c}        = [CONNECT_TEMPXC(c,1) CONNECT_TEMPXC(c,2) CONNECT_TEMPXC(c,3)];

     [N_localxc , dN_localxc] = test_new_bspline(Xi,deg,LOCC(cpid,1));
      
     for i = 1:3
         Nxc{c}(i) = N_localxc(CONNECT_TEMPXC(c,i));
         dNxc{c}(1,i) = dN_localxc(CONNECT_TEMPXC(c,i));
     end
    end
    
elseif choice ==2
%% Cubic Bspline
     deg = 3;
     Nelemsx = (NN(1)) - 4 + 2;
     XiCount = 1/le(1)+7; 
     Xi = zeros(XiCount,1);
     Xi(1:deg+1,1) = 2*le(1);
     for i=1:1/le(1)
     Xi(deg+1+i,1) =  Xi(deg+i,1)+le(1);
     end
     Xi(XiCount-deg+1:XiCount) = Xi(XiCount-deg);
     
    spElems         = zeros(spCount,1);
    CONNECT_TEMPX   = zeros(spCount,deg+1);
    
     for p = 1:spCount
     % compute vector store index elements 
     spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));

     % Knot vector index
     CONNECT_TEMPX(p,1) = floor(x_sp(p,1)/le(1)+1) - 2;
     CONNECT_TEMPX(p,2) = floor(x_sp(p,1)/le(1)+1) - 1;
     CONNECT_TEMPX(p,3) = floor(x_sp(p,1)/le(1)+1) - 0;
     CONNECT_TEMPX(p,4) = floor(x_sp(p,1)/le(1)+1) + 1;
     CONNECTX{p}        = [CONNECT_TEMPX(p,1) CONNECT_TEMPX(p,2) CONNECT_TEMPX(p,3) CONNECT_TEMPX(p,4)];
     [N_localx , dN_localx] = test_new_bspline(Xi,deg,x_sp(p,1));
     
     for i = 1:deg+1
         Nx{p}(i) = N_localx(CONNECT_TEMPX(p,i));
         dNx{p}(1,i) = dN_localx(CONNECT_TEMPX(p,i));
     end  
     end
     
     
     %% Compute bspline shape function from centroids to knot
    active_elements = unique(spElems);
    ElemsCount = length(active_elements);
    CONNECT_TEMPXC   = zeros(ElemsCount,3);

    for c = 1:ElemsCount
    cpid = active_elements(c);
    
    % Knot vector index for centroids 
     CONNECT_TEMPXC(c,1) = floor(LOCC(cpid,1)/le(1)+1) - 2;
     CONNECT_TEMPXC(c,2) = floor(LOCC(cpid,1)/le(1)+1) - 1;
     CONNECT_TEMPXC(c,3) = floor(LOCC(cpid,1)/le(1)+1) + 0;
     CONNECT_TEMPXC(c,4) = floor(LOCC(cpid,1)/le(1)+1) + 1;
     CONNECTXC{c}        = [CONNECT_TEMPXC(c,1) CONNECT_TEMPXC(c,2) CONNECT_TEMPXC(c,3) CONNECT_TEMPXC(c,4)];

     [N_localxc , dN_localxc] = test_new_bspline(Xi,deg,LOCC(cpid,1));
      
     for i = 1:4
         Nxc{c}(i) = N_localxc(CONNECT_TEMPXC(c,i));
         dNxc{c}(1,i) = dN_localxc(CONNECT_TEMPXC(c,i));
     end
    end
    
end
    
%% Compute the index of the node for the particles (copied)
    CONNECT_TEMP_G  = zeros(spCount,1);              % gauss number from left to right from bottom to top
    GAUSS           = 1 * ones(spCount,1);
   
     % index of the centroids to particles
     for p = 1:spCount
     % Centroid index for particle
%      CONNECT_TEMP_G(p,1) =  spElems(p)-(NN(1)-1)-1;
%      CONNECT_TEMP_G(p,2) =  spElems(p)-(NN(1)-1);
%      CONNECT_TEMP_G(p,3) =  spElems(p)-(NN(1)-1)+1;
%      CONNECT_TEMP_G(p,4) =  spElems(p)-1;
%      CONNECT_TEMP_G(p,5) =  spElems(p);
%      CONNECT_TEMP_G(p,6) =  spElems(p)+1;
%      CONNECT_TEMP_G(p,7) =  spElems(p)+(NN(1)-1)-1;
%      CONNECT_TEMP_G(p,8) =  spElems(p)+(NN(1)-1);
%      CONNECT_TEMP_G(p,9) =  spElems(p)+(NN(1)-1)+1;
%      CONNECT_G{p} = [CONNECT_TEMP_G(p,1) CONNECT_TEMP_G(p,2) CONNECT_TEMP_G(p,3) CONNECT_TEMP_G(p,4) CONNECT_TEMP_G(p,5) CONNECT_TEMP_G(p,6) CONNECT_TEMP_G(p,7) CONNECT_TEMP_G(p,8) CONNECT_TEMP_G(p,9)];
     
     CONNECT_TEMP_G(p,1)  = spElems(p);
     CONNECT_G{p} = spElems(p);
     end
     
     % Compute mspoints: index of particles in each element (active element)
 for c =1:cellCount
     id_p = find(spElems==c);
     mspoints{c}=id_p;
 end
 
% Least square matrix
Ac               = cell(cellCount,GaussCount);

% Function for mapping from particles to centroids using quadratic Bspline
Ng_local         = cell(spCount,1);                     % Value of quadratic bspline shape function
N_g              = cell(spCount,1);                     % Value of MLS shape function

for c = 1:cellCount
for g = 1:GaussCount
    Ac{c,g} = zeros(r,r);
end
end

for p = 1:spCount
    Ng_local{p} = zeros(GAUSS(p),GaussCount);
    N_g{p} = zeros(GAUSS(p),GaussCount);
end

% Compute the quadratic Bspline shape function for Gauss point
for p = 1:spCount
for c = 1:GAUSS(p)
    cpid = CONNECT_TEMP_G(p,c);

for g = 1:GaussCount
     % Compute the shape functions and gradient of the shape functions
        position_Gauss = [x_gauss(cpid,g)];
        dx = -x_sp(p,1)+position_Gauss;
        rx = abs(dx)/le(1);

        if (rx>1.5)
             Ng_local{p}(c,g)     = 0.0;
            elseif (rx<=0.5)
             Ng_local{p}(c,g)     = 3/4 - rx^2;
            else
             Ng_local{p}(c,g)     = 0.5 * (1.5 - rx)^2;
        end

 end
end
end

for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==2
         Pxy = [1 ; x_sp(p,1)];    
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,1)^2];   
     end
         
     PPxy = Pxy * Pxy';
     
     for c = 1:GAUSS(p)
         cpid = CONNECT_TEMP_G(p,c);
     for g = 1:GaussCount

         Ac{cpid,g} = Ac{cpid,g} + Ng_local{p}(c,g) * PPxy;
     end
    end
end
 
% Compute shape function and gradients
 for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==2
         Pxy = [1 ; x_sp(p,1)];   
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,1)^2];   
     end
     
     for c = 1:GAUSS(p)
         cpid = CONNECT_TEMP_G(p,c);
         
         for g = 1:GaussCount
         position_Gauss = [x_gauss(cpid,g)];

         if r==1
         pc = [1];
         elseif r==2
         pc = [1 ; position_Gauss(1)];   
         elseif r==3
         pc = [1 ; position_Gauss(1) ; position_Gauss(1)^2];   
         end
         
         if r==3
         
         rc = Ac{cpid,g} \ pc;
         N_g{p}(c,g) = rc' * Ng_local{p}(c,g) * Pxy;  
         
         elseif r==2
         rc = Ac{cpid,g}(1:2,1:2) \ pc(1:2);
         N_g{p}(c,g) = rc' * Ng_local{p}(c,g) * Pxy(1:2);  
         end
         end
     end
 end
 
%% Mapping from particle to nodes
% Interpolation from particle to grid task
% Node variables
nmass                = zeros(Nelemsx,1);                   % Nodal Mass
nmomentum            = zeros(Nelemsx,1);                   % Nodal Momentum
niforce              = zeros(Nelemsx,1);                   % Nodal Internal force
neforce              = zeros(Nelemsx,1);                   % Nodal External force
traction             = zeros(Nelemsx,1);                   % Nodal Traction
nforce               = zeros(Nelemsx,1);                   % Total force
nacc                 = zeros(Nelemsx,1);                   % Acc
nvelo                = zeros(Nelemsx,1);
% Nelemsx: number of basis function in x direction
% Nelemsy: number of basis function in y direction
% cell variables
p_sc                    = zeros(cellCount,GaussCount); 
pm_sc                   = zeros(cellCount,GaussCount); 
SSC                     = zeros(cellCount,GaussCount);
b_sc                    = zeros(cellCount,GaussCount); 

% Generate Gauss variables
for c = 1:cellCount
for g = 1:GaussCount
p_sc(c,g) = 0;
pm_sc(c,g) = 0;
SSC(c,g) = 0;
b_sc(c,g) = 0;
end
end

for p = 1 : spCount
%      % Build stress tensor
     SSP = [s_sp(p,1)];
         
     % Mapping from particle to centroids (copied)
     for c = 1:GAUSS(p)
     cpid = CONNECT_G{p}(c);
     if isempty(mspoints{cpid})==1
         continue
     end

    for g = 1:GaussCount
    position_Gauss = [x_gauss(cpid,g)];
    
    if position_Gauss(1)>2*le(1) && position_Gauss(1)<1+2*le(1)
     % Density
     p_sc(cpid,g)          = p_sc(cpid,g) + p_sp(p) * N_g{p}(c,g);
     % Momentum
     pm_sc(cpid,g)         = pm_sc(cpid,g) + p_sp(p) * v_ssp(p,1) * N_g{p}(c,g);
     % Body
     b_sc(cpid,g)          = b_sc(cpid,g) + p_sp(p) * b_sp(p,1) * N_g{p}(c,g);
     % Stress
     SSC(cpid,g)           = SSC(cpid,g) + SSP * N_g{p}(c,g);
    end
    
     end
    end
end
 
 %% Mapping from Gausses to knot
for c=1:ElemsCount
    cpid = active_elements(c);
     
    for g=1:GaussCount
        position_Gauss = [x_gauss(cpid,g)];
        
     [N_gaussx , dN_gaussx] = test_new_bspline(Xi,deg,position_Gauss(1));  

 for i = 1:deg+1
     xpid = CONNECTXC{c}(i); % global x coordinate of basis function
      
         Nxg = N_gaussx(xpid); dNxg = dN_gaussx(xpid);
         
         % Mass
         nmass(xpid)       = nmass(xpid) + p_sc(cpid,g) * Nxg * omega_gauss(cpid,g);
         
         % Momentum
         nmomentum(xpid)   = nmomentum(xpid) +  pm_sc(cpid,g) * Nxg * omega_gauss(cpid,g);

         % Internal force
         % Build stress tensor 
         niforce(xpid)     = niforce(xpid) - SSC(cpid,g) * dNxg * omega_gauss(c,g);
         
         % External force
         neforce(xpid)     = neforce(xpid) + b_sc(cpid,g) * Nxg * omega_gauss(cpid,g); 
         end
         end
end

% % Interpolation from particle to node
%  for p=1:spCount
%  % Build stress tensor
%  SSP = [s_sp(p,1)];
%  
%  for i = 1:deg+1
%      xpid = CONNECTX{p}(i); % global x coordinate of basis function        
%          % Mass
%          nmass(xpid)       = nmass(xpid)  + m_sp(p) * Nx{p}(i);
%          
%          % Momentum
%          nmomentum(xpid)    = nmomentum(xpid)  + m_sp(p) * v_ssp(p,1) * Nx{p}(i);
%      
%          % Internal force
%          niforce(xpid)      = niforce(xpid)  - V_sp(p) * SSP * dNx{p}(i);
%          
%          % External force
%          neforce(xpid)      = neforce(xpid)  + b_sp(p,1) * m_sp(p) * Nx{p}(i);         
%      end
%  end

Gauss_position = [];
data_density = [];
for c = 1:length(active_elements)  
    cpid = active_elements(c);
    Gauss_position = [Gauss_position x_gauss(cpid,1:2)];
    data_density = [data_density p_sc(cpid,1:2)*(le(1)/GaussCount)];
end

% Consistent Mass
 [N_consistent , ~] = test_new_bspline(Xi,deg,Gauss_position');
 
 MASS = bsxfun(@times, data_density', N_consistent);
 MASS = MASS'*N_consistent;
 MASS = sparse(MASS);
 
 MASS_diagonal = trace(MASS)/length(nmomentum);
 boundary_knot = [1 length(nmomentum)];
 MASS(:,boundary_knot) = 0;
 MASS(boundary_knot,:) = 0;
 MASS(boundary_knot,boundary_knot) = MASS_diagonal*speye(length(boundary_knot));

 if lump==1
 % Lump mass
 MASS = zeros(Nelemsx,Nelemsx);
 for i = 1:Nelemsx
     MASS(i,i) = nmass(i);
 end
 end
%% Update momentum
% Update force and momentum
for i = 1:Nelemsx
        nforce(i)     	= niforce(i) + neforce(i) + traction(i);
%         nacc(i)         = nforce(i)/nmass(i);
        nvelo(i)        = nmomentum(i)/nmass(i);
end

nacc = MASS\nforce;

for i = 1:Nelemsx
%         nmomentum(i)      = nmomentum(i) + nforce(i)*dt;
        nvelo(i)      = nvelo(i) + nacc(i)*dt;
        % Boundary condition
        if i == 1 || i == Nelemsx
%             nforce(i)      = 0;
%             nmomentum(i)   = 0;
            nvelo(i)       = 0;
            nacc(i)       = 0;
        end
end

%% Update solid particle velocity and position

for p = 1:spCount
    for i = 1:deg+1
     xpid = CONNECTX{p}(i); % global x coordinate of basis function         
         if nmass(i) ==0
             continue
         end
         
%          v_ssp(p,1)                      = v_ssp(p,1) + dt * Nx{p}(i) * nforce(xpid)/nmass(xpid);
%          x_sp(p,1)                       = x_sp(p,1) + dt * Nx{p}(i) * nmomentum(xpid)/nmass(xpid);
         v_ssp(p,1)                      = v_ssp(p,1) + dt * Nx{p}(i) * nacc(xpid);
         x_sp(p,1)                       = x_sp(p,1) + dt * Nx{p}(i) * nvelo(xpid);
         d_sp(p,1)                       = x_sp(p,1) - x_spo(p,1);
    end
end
     
% %% Mapping nodal velocity back to node
% % Node variables
%  nvelo                = zeros(Nelemsx,1);                   % Nodal Velocity
%  nmomentum            = zeros(Nelemsx,1);                   % Nodal Momentum
%  
%  for i = 1:Nelemsx
%         nmomentum(i)         = 0;
%         nvelo(i)             = 0;  
%  end
% 
%  % Interpolation momentum from particle to node
%  for p=1:spCount 
%  for i = 1:deg+1
%      xpid = CONNECTX{p}(i); % global x coordinate of basis function        
%          % Momentum
%          nmomentum(xpid)   = nmomentum(xpid) + m_sp(p) * v_ssp(p,1) * Nx{p}(i);
%  end
%  end
% 
%  nvelo = MASS\nmomentum;
%  
%  % Velocity
%  for i = 1:Nelemsx
%      if  nmass(i) == 0
%          continue
%      end
%      
% %         nvelo(i)     	= nmomentum(i) / nmass(i);
%         
%         % Boundary condition
%         if i == 1 || i == Nelemsx
%             nvelo(i)      = 0;
%         end
%  end

% if choice ==1
% %% Quadratic Bspline
%      deg = 2;
%      Nelemsx = (NN(1)) - 4 + 1;
%      XiCount = 1/le(1)+5; 
%      Xi = zeros(XiCount,1);
%      Xi(1:deg+1,1) = 2*le(1);
%      for i=1:1/le(1)
%      Xi(deg+1+i,1) =  Xi(deg+i,1)+le(1);
%      end
%      Xi(XiCount-1:XiCount) = Xi(XiCount-2);
%      
%     spElems         = zeros(spCount,1);
%     CONNECT_TEMPX   = zeros(spCount,deg+1);
%      for p = 1:spCount
%      % compute vector store index elements 
%      spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));
% 
%      % Knot vector index
%      CONNECT_TEMPX(p,1) = floor(x_sp(p,1)/le(1)+1) - 2;
%      CONNECT_TEMPX(p,2) = floor(x_sp(p,1)/le(1)+1) - 1;
%      CONNECT_TEMPX(p,3) = floor(x_sp(p,1)/le(1)+1) - 0;
%      CONNECTX{p}        = [CONNECT_TEMPX(p,1) CONNECT_TEMPX(p,2) CONNECT_TEMPX(p,3)];
%      [N_localx , dN_localx] = test_new_bspline(Xi,deg,x_sp(p,1));
%      
%      for i = 1:deg+1
%          Nx{p}(i) = N_localx(CONNECT_TEMPX(p,i));
%          dNx{p}(1,i) = dN_localx(CONNECT_TEMPX(p,i));
%      end  
%      end
% 
% elseif choice ==2
% %% Cubic Bspline
%      deg = 3;
%      Nelemsx = (NN(1)) - 4 + 6;
%      XiCount = 1/le(1)+7; 
%      Xi = zeros(XiCount,1);
%      Xi(1:deg+1,1) = 2*le(1);
%      for i=1:1/le(1)
%      Xi(deg+1+i,1) =  Xi(deg+i,1)+le(1);
%      end
%      Xi(XiCount-deg+1:XiCount) = Xi(XiCount-deg);
%      
%     spElems         = zeros(spCount,1);
%     CONNECT_TEMPX   = zeros(spCount,deg+1);
%     
%      for p = 1:spCount
%      % compute vector store index elements 
%      spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));
% 
%      % Knot vector index
%      CONNECT_TEMPX(p,1) = floor(x_sp(p,1)/le(1)+1) - 2;
%      CONNECT_TEMPX(p,2) = floor(x_sp(p,1)/le(1)+1) - 1;
%      CONNECT_TEMPX(p,3) = floor(x_sp(p,1)/le(1)+1) - 0;
%      CONNECT_TEMPX(p,4) = floor(x_sp(p,1)/le(1)+1) + 1;
%      CONNECTX{p}        = [CONNECT_TEMPX(p,1) CONNECT_TEMPX(p,2) CONNECT_TEMPX(p,3) CONNECT_TEMPX(p,4)];
%      [N_localx , dN_localx] = test_new_bspline(Xi,deg,x_sp(p,1));
%      
%      for i = 1:deg+1
%          Nx{p}(i) = N_localx(CONNECT_TEMPX(p,i));
%          dNx{p}(1,i) = dN_localx(CONNECT_TEMPX(p,i));
%      end  
%      end
%      
% end

%% Update effective stress
L_sp = cell(spCount,1);

for sp = 1:spCount
    L_sp{sp} = 0;
end


for p = 1:spCount
    for i = 1:deg+1
     xpid = CONNECTX{p}(i); % global x coordinate of basis function      
        L_sp{p}   = L_sp{p} + nvelo(xpid) * dNx{p}(i);
    end        
        F_sp{p}(1,1) = (1+L_sp{p}*dt)*F_sp{p}(1,1);                           
        J = det(F_sp{p}(1,1));
        V_sp(p)=V_spo(p)*J;   

        switch CModel
            case 'Neo_Hookean_Elastic'
                [s_sp(p,:)]=Neo_Hookean_elastic(CModel_parameter,F_sp{p},J);
            case 'Linear_Elastic'
                [s_sp(p,:)]=Linear_elastic(CModel_parameter,dESP,s_sp(p,:));             
            case 'Water'
                [s_sp(p,:)]=Water(CModel_parameter,(L_sp{p} + L_sp{p}')/2,J);
        end
        p_sp(p) = m_sp(p)/V_sp(p);
end      
