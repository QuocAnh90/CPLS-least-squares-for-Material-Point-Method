function[v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = BSMPM_solvertest(CModel,CModel_parameter,...
    nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
    nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt)

%% Store particles into cell
choice =1;

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

% Interpolation from particle to node
 for p=1:spCount
 % Build stress tensor
 SSP = [s_sp(p,1)];
 
 for i = 1:deg+1
     xpid = CONNECTX{p}(i); % global x coordinate of basis function        
         % Mass
         nmass(xpid)       = nmass(xpid)  + m_sp(p) * Nx{p}(i);
         
         % Momentum
         nmomentum(xpid)    = nmomentum(xpid)  + m_sp(p) * v_ssp(p,1) * Nx{p}(i);
     
         % Internal force
         niforce(xpid)      = niforce(xpid)  - V_sp(p) * SSP * dNx{p}(i);
         
         % External force
         neforce(xpid)      = neforce(xpid)  + b_sp(p,1) * m_sp(p) * Nx{p}(i);         
     end
 end

% Consistent Mass
 [N_consistent , ~] = test_new_bspline(Xi,deg,x_sp(:,1));
 
 MASS = bsxfun(@times, m_sp, N_consistent);
 MASS = MASS'*N_consistent;
 MASS = sparse(MASS);
 
 MASS_diagonal = trace(MASS)/length(nmomentum);
 boundary_knot = [1 length(nmomentum)];
 MASS(:,boundary_knot) = 0;
 MASS(boundary_knot,:) = 0;
 MASS(boundary_knot,boundary_knot) = MASS_diagonal*speye(length(boundary_knot));

 % Lump mass
 MASS = zeros(Nelemsx,Nelemsx);
 for i = 1:Nelemsx
     MASS(i,i) = nmass(i);
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
