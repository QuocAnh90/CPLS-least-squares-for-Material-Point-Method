function[v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MLS_BSMPM_solver(CModel,CModel_parameter,...
    nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOCC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
    F_sp,V_spo,dt,GaussCount,x_gauss,y_gauss,omega_gauss)

%% Store particles into cell
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

%% Parameters
spElems         = zeros(spCount,1);              % index of elements where stores particles
CONNECT_TEMPX   = zeros(spCount,3);
CONNECT_TEMPY   = zeros(spCount,3); 
CONNECT_TEMP_G  = zeros(spCount,9);              % gauss number from left to right from bottom to top
CONNECT_G       = cell(spCount,1);
GAUSS           = 9 * ones(spCount,1);
NODES           = 4 * ones(spCount,1);

% Order of MLS
r = 3;
r_c = r * ones(cellCount,1);
r_n = r * ones(nodeCount,1);

 % Centroid index for particle
 for p = 1:spCount
 spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));
 
  % Node index for particle using linear basis
 CONNECT_TEMP(p,1) = spElems(p) + ceil(spElems(p)/(NN(1)-1)) - 1;
 CONNECT_TEMP(p,2) = CONNECT_TEMP(p,1)+1; 
 CONNECT_TEMP(p,3) = CONNECT_TEMP(p,2)+NN(1); 
 CONNECT_TEMP(p,4) = CONNECT_TEMP(p,1)+NN(1);
 CONNECT{p}        = [CONNECT_TEMP(p,1) CONNECT_TEMP(p,2) CONNECT_TEMP(p,3) CONNECT_TEMP(p,4)];
 
 end
 
  % Compute mspoints: index of particles in each element (active element)
 for c =1:cellCount
     id_p = find(spElems==c);
     mspoints{c}=id_p;
 end
 
%% Knot vector reconstruction (copied)
     deg = 2;
     Nelemsx = (NN(1)) - 4 + 1;
     Nelemsy = 3;
%      Nelemsy = (NN(2)) - 4 + 1;
     
     XiCount = 1/le(1)+5; 
     Xi = zeros(XiCount,1);
     Xi(1:3,1) = 2*le(1);
     for i=1:1/le(1)
     Xi(3+i,1) =  Xi(2+i,1)+le(1);
     end
     Xi(XiCount-1:XiCount) = Xi(XiCount-2);
     
     YiCount = 3;
     Yi = zeros(YiCount,1);
     Yi(1:3,1) = 1*le(2);
     Yi(4:6,1) = 2*le(2);

%      YiCount = XiCount;
%      Yi = Xi;
     
 %% Compute bspline shape function from particle to knot (copied)
      for p = 1:spCount
     % compute vector store index elements 
     spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));

     % Knot vector index
     CONNECT_TEMPX(p,1) = floor(x_sp(p,1)/le(1)+1) - 2;
     CONNECT_TEMPX(p,2) = floor(x_sp(p,1)/le(1)+1) - 1;
     CONNECT_TEMPX(p,3) = floor(x_sp(p,1)/le(1)+1) - 0;
     CONNECTX{p}        = [CONNECT_TEMPX(p,1) CONNECT_TEMPX(p,2) CONNECT_TEMPX(p,3)];

     CONNECT_TEMPY(p,1) = floor(x_sp(p,2)/le(2)+1) - 1;
     CONNECT_TEMPY(p,2) = floor(x_sp(p,2)/le(2)+1) - 0;
     CONNECT_TEMPY(p,3) = floor(x_sp(p,2)/le(2)+1) + 1;
     CONNECTY{p}        = [CONNECT_TEMPY(p,1) CONNECT_TEMPY(p,2) CONNECT_TEMPY(p,3)];
 
     [N_localx , dN_localx] = test_new_bspline(Xi,deg,x_sp(p,1));
     
     for i = 1:3
         Nx{p}(i) = N_localx(CONNECT_TEMPX(p,i));
         dNx{p}(1,i) = dN_localx(CONNECT_TEMPX(p,i));
     end
     
     [N_localy , dN_localy] = test_new_bspline(Yi,deg,x_sp(p,2));
     
     for i = 1:3
         Ny{p}(i) = N_localy(CONNECT_TEMPY(p,i));
         dNy{p}(1,i) = dN_localy(CONNECT_TEMPY(p,i));
     end
     end
 
 %% Compute bspline shape function from centroids to knot
    active_elements = unique(spElems);
    ElemsCount = length(active_elements);
    CONNECT_TEMPXC   = zeros(ElemsCount,3);
    CONNECT_TEMPYC   = zeros(ElemsCount,3); 

    for c = 1:ElemsCount
    cpid = active_elements(c);
    
    % Knot vector index for centroids 
     CONNECT_TEMPXC(c,1) = floor(LOCC(cpid,1)/le(1)+1) - 2;
     CONNECT_TEMPXC(c,2) = floor(LOCC(cpid,1)/le(1)+1) - 1;
     CONNECT_TEMPXC(c,3) = floor(LOCC(cpid,1)/le(1)+1) + 0;
     CONNECTXC{c}        = [CONNECT_TEMPXC(c,1) CONNECT_TEMPXC(c,2) CONNECT_TEMPXC(c,3)];

     CONNECT_TEMPYC(c,1) = floor(LOCC(cpid,2)/le(2)+1) - 1;
     CONNECT_TEMPYC(c,2) = floor(LOCC(cpid,2)/le(2)+1) - 0;
     CONNECT_TEMPYC(c,3) = floor(LOCC(cpid,2)/le(2)+1) + 1;
     CONNECTYC{c}        = [CONNECT_TEMPYC(c,1) CONNECT_TEMPYC(c,2) CONNECT_TEMPYC(c,3)];
 
%      [N_localxc , dN_localxc] = test_new_bspline(Xi,deg,LOCC(cpid,1));
%      
%      for i = 1:3
%          Nxc{c}(i) = N_localxc(CONNECT_TEMPXC(c,i));
%          dNxc{c}(1,i) = dN_localxc(CONNECT_TEMPXC(c,i));
%      end
%      
%      [N_localyc , dN_localyc] = test_new_bspline(Yi,deg,LOCC(cpid,2));
%      
%      for i = 1:3
%          Nyc{c}(i) = N_localyc(CONNECT_TEMPYC(c,i));
%          dNyc{c}(1,i) = dN_localyc(CONNECT_TEMPYC(c,i));
%      end   
    end
    
 %% Compute the index of the node for the particles (copied)
     % index of the centroids to particles
     for p = 1:spCount
     % Centroid index for particle
     CONNECT_TEMP_G(p,1) =  spElems(p)-(NN(1)-1)-1;
     CONNECT_TEMP_G(p,2) =  spElems(p)-(NN(1)-1);
     CONNECT_TEMP_G(p,3) =  spElems(p)-(NN(1)-1)+1;
     CONNECT_TEMP_G(p,4) =  spElems(p)-1;
     CONNECT_TEMP_G(p,5) =  spElems(p);
     CONNECT_TEMP_G(p,6) =  spElems(p)+1;
     CONNECT_TEMP_G(p,7) =  spElems(p)+(NN(1)-1)-1;
     CONNECT_TEMP_G(p,8) =  spElems(p)+(NN(1)-1);
     CONNECT_TEMP_G(p,9) =  spElems(p)+(NN(1)-1)+1;
     CONNECT_G{p} = [CONNECT_TEMP_G(p,1) CONNECT_TEMP_G(p,2) CONNECT_TEMP_G(p,3) CONNECT_TEMP_G(p,4) CONNECT_TEMP_G(p,5) CONNECT_TEMP_G(p,6) CONNECT_TEMP_G(p,7) CONNECT_TEMP_G(p,8) CONNECT_TEMP_G(p,9)];
     
%     CONNECT_TEMP_G  = zeros(spCount,1);              % gauss number from left to right from bottom to top
%     GAUSS           = 1 * ones(spCount,1);
%      CONNECT_TEMP_G(p,1)  = spElems(p);
%      CONNECT_G{p} = spElems(p);
     end
 
 %% Quadratic-Bspline for mapping from particle to centroids using MLS
   N_g = Quadratic4Gauss_MLS(cellCount,spCount,GAUSS,CONNECT_TEMP_G,x_sp,GaussCount,x_gauss,y_gauss,le,r,r_c);
   
%% Mapping from particle to centroids
% cell variables
p_sc                    = cell(cellCount,GaussCount); 
pm_sc                   = cell(cellCount,GaussCount); 
SSC                     = cell(cellCount,GaussCount);
b_sc                    = cell(cellCount,GaussCount); 

% Generate Gauss variables
for c = 1:cellCount
for g = 1:GaussCount
p_sc{c,g} = 0;
pm_sc{c,g} = zeros(1,2);
SSC{c,g} = zeros(2,2);
b_sc{c,g} = zeros(1,2);
end

end

% Node variables
nmass                   = cell(Nelemsx,Nelemsy);                   % Nodal Mass
nmomentum               = cell(Nelemsx,Nelemsy);                   % Nodal Momentum
niforce                 = cell(Nelemsx,Nelemsy);                   % Nodal Internal force
neforce                 = cell(Nelemsx,Nelemsy);                   % Nodal External force
traction                = cell(Nelemsx,Nelemsy);                   % Nodal Traction
nforce                  = cell(Nelemsx,Nelemsy);                   % Total force

% Generate the nodal data
for i = 1:Nelemsx
    for j = 1:Nelemsy
        nmass{i,j}                = 0;                   
        nmomentum{i,j}            = zeros(1,2);
        niforce{i,j}              = zeros(1,2);  
        neforce{i,j}              = zeros(1,2);    
        traction{i,j}             = zeros(1,2);      
    end
end

for p = 1 : spCount
%      % Build stress tensor
     SSP = [s_sp(p,1) s_sp(p,3);s_sp(p,3) s_sp(p,2)];
         
     % Mapping from particle to centroids (copied)
     for c = 1:GAUSS(p)
     cpid = CONNECT_G{p}(c);
     if isempty(mspoints{cpid})==1
         continue
     end

    for g = 1:GaussCount

     % Density
     p_sc{cpid,g}          = p_sc{cpid,g} + p_sp(p) * N_g{p}(c,g);
     % Momentum
     pm_sc{cpid,g}         = pm_sc{cpid,g} + p_sp(p) * v_ssp(p,:) * N_g{p}(c,g);
     % Body
     b_sc{cpid,g}          = b_sc{cpid,g} + p_sp(p) * b_sp(p,:) * N_g{p}(c,g);
     % Stress
     SSC{cpid,g}           = SSC{cpid,g} + SSP * N_g{p}(c,g);
     end
    end
 end
 
 %% Mapping from Gausses to knot
for c=1:ElemsCount
    cpid = active_elements(c);
     
    for g=1:GaussCount
        position_Gauss = [x_gauss(cpid,g);y_gauss(cpid,g)];
        
     [N_gaussx , dN_gaussx] = test_new_bspline(Xi,deg,position_Gauss(1));  
     [N_gaussy , dN_gaussy] = test_new_bspline(Yi,deg,position_Gauss(2));

 for i = 1:3
     xpid = CONNECTXC{c}(i); % global x coordinate of basis function
     for j = 1:3
         ypid = CONNECTYC{c}(j); % global y coordinate of basis function
         
         Nxg = N_gaussx(xpid); dNxg = dN_gaussx(xpid);
         Nyg = N_gaussy(ypid); dNyg = dN_gaussy(ypid);
         
         % Mass
         nmass{xpid,ypid}       = nmass{xpid,ypid} + p_sc{cpid,g} * Nxg * Nyg * omega_gauss(cpid,g);
         
         % Momentum
         nmomentum{xpid,ypid}   = nmomentum{xpid,ypid} +  pm_sc{cpid,g} * Nxg * Nyg * omega_gauss(cpid,g);

         % Internal force
         % Build stress tensor 
         niforce{xpid,ypid}     = niforce{xpid,ypid} - (SSC{cpid,g} * [dNxg*Nyg;Nxg*dNyg])'*omega_gauss(c,g);
         
         % External force
         neforce{xpid,ypid}     = neforce{xpid,ypid} + b_sc{cpid,g} * Nxg * Nyg * omega_gauss(cpid,g); 
     end
 end
    end
 end

%% Update momentum (copied)
% Update force and momentum
for i = 1:Nelemsx
    for j = 1:Nelemsy
        nforce{i,j}     	= niforce{i,j} + neforce{i,j} + traction{i,j};
    end
end

for i = 1:Nelemsx
    for j = 1:Nelemsy
        nmomentum{i,j}      = nmomentum{i,j} + nforce{i,j}*dt;
        
        % Boundary condition
        if i == 1 || i == Nelemsx
            nforce{i,j}(1)      = 0;
            nmomentum{i,j}(1)   = 0;
        end
        
        if j == 1 || j == Nelemsy
            nforce{i,j}(2)      = 0;
            nmomentum{i,j}(2)   = 0;
        end
    end
end

%% Update solid particle velocity and position (copied)

for p = 1:spCount
    for i = 1:3
     xpid = CONNECTX{p}(i); % global x coordinate of basis function
    for j = 1:3
         ypid = CONNECTY{p}(j); % global y coordinate of basis function
         
         if nmass{i,j} ==0
             continue
         end
         
         v_ssp(p,:)                      = v_ssp(p,:) + dt * Nx{p}(i) * Ny{p}(j) * nforce{xpid,ypid}/nmass{xpid,ypid};
         x_sp(p,:)                       = x_sp(p,:) + dt * Nx{p}(i) * Ny{p}(j) * nmomentum{xpid,ypid}/nmass{xpid,ypid};
         d_sp(p,:)                       = x_sp(p,:) - x_spo(p,:);
    end
    end
end

%  %% Compute bspline shape function from particle to knot (copied)
%       for p = 1:spCount
%      % compute vector store index elements 
%      spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));
% 
%      % Knot vector index
%      CONNECT_TEMPX(p,1) = floor(x_sp(p,1)/le(1)+1) - 2;
%      CONNECT_TEMPX(p,2) = floor(x_sp(p,1)/le(1)+1) - 1;
%      CONNECT_TEMPX(p,3) = floor(x_sp(p,1)/le(1)+1) - 0;
%      CONNECTX{p}        = [CONNECT_TEMPX(p,1) CONNECT_TEMPX(p,2) CONNECT_TEMPX(p,3)];
% 
%      CONNECT_TEMPY(p,1) = floor(x_sp(p,2)/le(2)+1) - 1;
%      CONNECT_TEMPY(p,2) = floor(x_sp(p,2)/le(2)+1) - 0;
%      CONNECT_TEMPY(p,3) = floor(x_sp(p,2)/le(2)+1) + 1;
%      CONNECTY{p}        = [CONNECT_TEMPY(p,1) CONNECT_TEMPY(p,2) CONNECT_TEMPY(p,3)];
%  
%      [N_localx , dN_localx] = test_new_bspline(Xi,deg,x_sp(p,1));
%      
%      for i = 1:3
%          Nx{p}(i) = N_localx(CONNECT_TEMPX(p,i));
%          dNx{p}(1,i) = dN_localx(CONNECT_TEMPX(p,i));
%      end
%      
%      [N_localy , dN_localy] = test_new_bspline(Yi,deg,x_sp(p,2));
%      
%      for i = 1:3
%          Ny{p}(i) = N_localy(CONNECT_TEMPY(p,i));
%          dNy{p}(1,i) = dN_localy(CONNECT_TEMPY(p,i));
%      end
%       end
%      
%        %% Compute the index of the node for the particles (copied)
%      % index of the centroids to particles
%      for p = 1:spCount
%      % Centroid index for particle
%      CONNECT_TEMP_G(p,1) =  spElems(p)-(NN(1)-1)-1;
%      CONNECT_TEMP_G(p,2) =  spElems(p)-(NN(1)-1);
%      CONNECT_TEMP_G(p,3) =  spElems(p)-(NN(1)-1)+1;
%      CONNECT_TEMP_G(p,4) =  spElems(p)-1;
%      CONNECT_TEMP_G(p,5) =  spElems(p);
%      CONNECT_TEMP_G(p,6) =  spElems(p)+1;
%      CONNECT_TEMP_G(p spElems(p);
%      CONNECT_G{p} = spElems(,7) =  spElems(p)+(NN(1)-1)-1;
%      CONNECT_TEMP_G(p,8) =  spElems(p)+(NN(1)-1);
%      CONNECT_TEMP_G(p,9) =  spElems(p)+(NN(1)-1)+1;
%      CONNECT_G{p} = [CONNECT_TEMP_G(p,1) CONNECT_TEMP_G(p,2) CONNECT_TEMP_G(p,3) CONNECT_TEMP_G(p,4) CONNECT_TEMP_G(p,5) CONNECT_TEMP_G(p,6) CONNECT_TEMP_G(p,7) CONNECT_TEMP_G(p,8) CONNECT_TEMP_G(p,9)];
%   
% %     CONNECT_TEMP_G  = zeros(spCount,1);              % gauss number from left to right from bottom to top
% %     GAUSS           = 1 * ones(spCount,1);
% %    
% %      CONNECT_TEMP_G(p,1)  = spElems(p);
% %      CONNECT_G{p} = spElems(p);
%      end
%      
%        %% Quadratic-Bspline for mapping from particle to centroids using MLS
%    N_g = Quadratic4Gauss_MLS(cellCount,spCount,GAUSS,CONNECT_TEMP_G,x_sp,GaussCount,x_gauss,y_gauss,le,r,r_c);
  
   
%% Mapping nodal velocity back to knot
pm_sc                   = cell(cellCount,GaussCount); 

% Generate Gauss variables
for c = 1:cellCount
for g = 1:GaussCount
pm_sc{c,g} = zeros(1,2);
end

end

% Node variables
 nvelo                = cell(Nelemsx,Nelemsy);                   % Nodal Velocity
 nmomentum            = cell(Nelemsx,Nelemsy);                   % Nodal Momentum

% Generate the nodal data
for i = 1:Nelemsx
    for j = 1:Nelemsy
        nmomentum{i,j}            = zeros(1,2);
        nvelo{i,j}                = zeros(1,2);
    end
end

for p = 1 : spCount
     % Mapping from particle to centroids (copied)
     for c = 1:GAUSS(p)
     cpid = CONNECT_G{p}(c);
     if isempty(mspoints{cpid})==1
         continue
     end

    for g = 1:GaussCount
     % Momentum
     pm_sc{cpid,g}         = pm_sc{cpid,g} + p_sp(p) * v_ssp(p,:) * N_g{p}(c,g);
     end
    end
 end
 
 %% Mapping from Gausses to knot
for c=1:ElemsCount
    cpid = active_elements(c);
     
    for g=1:GaussCount
        position_Gauss = [x_gauss(cpid,g);y_gauss(cpid,g)];
        
     [N_gaussx , dN_gaussx] = test_new_bspline(Xi,deg,position_Gauss(1));  
     [N_gaussy , dN_gaussy] = test_new_bspline(Yi,deg,position_Gauss(2));

 for i = 1:3
     xpid = CONNECTXC{c}(i); % global x coordinate of basis function
     for j = 1:3
         ypid = CONNECTYC{c}(j); % global y coordinate of basis function
         
         Nxg = N_gaussx(xpid); dNxg = dN_gaussx(xpid);
         Nyg = N_gaussy(ypid); dNyg = dN_gaussy(ypid);

         % Momentum
         nmomentum{xpid,ypid}   = nmomentum{xpid,ypid} +  pm_sc{cpid,g} * Nxg * Nyg * omega_gauss(cpid,g);
     end
 end
    end
 end

 % Velocity
 for i = 1:Nelemsx
    for j = 1:Nelemsy
        if nmass{i,j} ==0
            continue
        end
        
        nvelo{i,j}     	= nmomentum{i,j} / nmass{i,j};
        
        % Boundary condition
        if i == 1 || i == Nelemsx
            nvelo{i,j}(1)      = 0;
        end
        
        if j == 1 || j == Nelemsy
            nvelo{i,j}(2)      = 0;
        end
    end
 end
 
% figure
% plot(LOCC(:,1),b_sc(:,1)* le(1)*le(2))
% figure
% plot(x_sp(:,1),b_sp(:,1))

%% Update effective stress (copied)
L_sp = cell(spCount,1);

for sp = 1:spCount
    L_sp{sp} = zeros(2,2);
end


for p = 1:spCount
    for i = 1:3
     xpid = CONNECTX{p}(i); % global x coordinate of basis function
    for j = 1:3
         ypid = CONNECTY{p}(j); % global y coordinate of basis function
         
         L_sp{p}   = L_sp{p} + nvelo{xpid,ypid}' * [dNx{p}(i)*Ny{p}(j) Nx{p}(i)*dNy{p}(j)];
    end
    end
    dESP = (L_sp{p} + L_sp{p}')/2*dt; 
        
        F_sp{p} = (eye(2,2)+L_sp{p}*dt)*F_sp{p};                           
        J = det(F_sp{p});
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