function [k,node_boundary,color_cell,color_node,N_MLS,...
    CONNECT_MLS]...
    = Compute_Interpolator_IMLS_CPDIv2(x_data,...
    nodeCount,spCount,cellCount,x_sp,le,NN,LOC,r1_sp,r2_sp)

%% CPDI inputs
%  CONNECT_TEMP_G = zeros(spCount,36);
%  CONNECT_sp     = zeros(spCount,4);
 CONNECT_sp     = zeros(spCount,16);

%  CONNECT_sp_G   = zeros(spCount,9);
 spElems_corner = zeros(spCount,4);
 spElems        = zeros(spCount,1);
 npoints        = cell(nodeCount,1);
%  cpoints        = cell(cellCount,1);
%  GAUSS          = 9 * ones(spCount,1);
 msdata         = cell(cellCount,1);
%  cElems         = zeros(cellCount,1);
 CONNECT_MLS    = cell(spCount,5);
%  CONNECT_MLS_G  = cell(spCount,5); 
 
 color_cell     = zeros(cellCount,1);

 % Compute corner positions
x_corner1 = x_sp - r1_sp - r2_sp;
x_corner2 = x_sp + r1_sp - r2_sp;
x_corner3 = x_sp + r1_sp + r2_sp;
x_corner4 = x_sp - r1_sp + r2_sp;
 
 %% Build connectivity
 
 for sp = 1:spCount
 spElems(sp) = ceil(x_sp(sp,1)/le(1))+(NN(1)-1)*(fix(x_sp(sp,2)/le(2)));   
 % compute vector store elements index
%  msdata{spElems(sp)} = [msdata{spElems(sp)} [sp,1]];
%  % compute vector store data index
%  
 spElems_corner(sp,1) = ceil(x_corner1(sp,1)/le(1))+(NN(1)-1)*(fix(x_corner1(sp,2)/le(2)));                        
%  msdata{spElems_corner(sp,1)} = [msdata{spElems_corner(sp,1)} [sp,2]];
 spElems_corner(sp,2) = ceil(x_corner2(sp,1)/le(1))+(NN(1)-1)*(fix(x_corner2(sp,2)/le(2)));                        
%  msdata{spElems_corner(sp,2)} = [msdata{spElems_corner(sp,2)} [sp,3]];
 spElems_corner(sp,3) = ceil(x_corner3(sp,1)/le(1))+(NN(1)-1)*(fix(x_corner3(sp,2)/le(2)));
%  msdata{spElems_corner(sp,3)} = [msdata{spElems_corner(sp,3)} [sp,4]];
 spElems_corner(sp,4) = ceil(x_corner4(sp,1)/le(1))+(NN(1)-1)*(fix(x_corner4(sp,2)/le(2)));
%  msdata{spElems_corner(sp,4)} = [msdata{spElems_corner(sp,4)} [sp,5]];

 % Particle-node interaction
%  CONNECT_sp(sp,1) = spElems(sp) + ceil(spElems(sp)/(NN(1)-1)) - 1;
%  npoints{CONNECT_sp(sp,1)}   = [npoints{CONNECT_sp(sp,1)} [sp;1]];
%  CONNECT_sp(sp,2) = CONNECT_sp(sp,1)+1;
%  npoints{CONNECT_sp(sp,2)}   = [npoints{CONNECT_sp(sp,2)} [sp;1]];
%  CONNECT_sp(sp,3) = CONNECT_sp(sp,2)+NN(1);
%  npoints{CONNECT_sp(sp,3)}   = [npoints{CONNECT_sp(sp,3)} [sp;1]];
%  CONNECT_sp(sp,4) = CONNECT_sp(sp,1)+NN(1);
%  npoints{CONNECT_sp(sp,4)}   = [npoints{CONNECT_sp(sp,4)} [sp;1]];
 
 CONNECT_sp(sp,6) = spElems(sp) + ceil(spElems(sp)/(NN(1)-1)) - 1;
 npoints{CONNECT_sp(sp,6)}   = [npoints{CONNECT_sp(sp,6)} [sp;1]];
 CONNECT_sp(sp,7) = CONNECT_sp(sp,6)+1;
 npoints{CONNECT_sp(sp,7)}   = [npoints{CONNECT_sp(sp,7)} [sp;1]];
 CONNECT_sp(sp,11) = CONNECT_sp(sp,7)+NN(1);
 npoints{CONNECT_sp(sp,11)}   = [npoints{CONNECT_sp(sp,11)} [sp;1]];
 CONNECT_sp(sp,10) = CONNECT_sp(sp,6)+NN(1);
 npoints{CONNECT_sp(sp,10)}   = [npoints{CONNECT_sp(sp,10)} [sp;1]];
 
 CONNECT_sp(sp,1) = CONNECT_sp(sp,6) - NN(1) - 1;                 
 npoints{CONNECT_sp(sp,1)} = [npoints{CONNECT_sp(sp,1)} [sp;1]];
 CONNECT_sp(sp,2) = CONNECT_sp(sp,1) + 1;                         
 npoints{CONNECT_sp(sp,2)} = [npoints{CONNECT_sp(sp,2)} [sp;1]];
 CONNECT_sp(sp,3) = CONNECT_sp(sp,2) + 1;                         
 npoints{CONNECT_sp(sp,3)} = [npoints{CONNECT_sp(sp,3)} [sp;1]];
 CONNECT_sp(sp,4) = CONNECT_sp(sp,3) + 1;                         
 npoints{CONNECT_sp(sp,4)} = [npoints{CONNECT_sp(sp,4)} [sp;1]];
 
 CONNECT_sp(sp,5) = CONNECT_sp(sp,6) - 1;                         
 npoints{CONNECT_sp(sp,5)} = [npoints{CONNECT_sp(sp,5)} [sp;1]];
 CONNECT_sp(sp,8) = CONNECT_sp(sp,7) + 1;                         
 npoints{CONNECT_sp(sp,8)} = [npoints{CONNECT_sp(sp,8)} [sp;1]];

CONNECT_sp(sp,9) = CONNECT_sp(sp,5) + NN(1);                     
npoints{CONNECT_sp(sp,9)} = [npoints{CONNECT_sp(sp,9)} [sp;1]];
 CONNECT_sp(sp,12)= CONNECT_sp(sp,11) + 1;                        
npoints{CONNECT_sp(sp,12)} = [npoints{CONNECT_sp(sp,12)} [sp;1]];

 CONNECT_sp(sp,13)= CONNECT_sp(sp,9) + NN(1);                     
npoints{CONNECT_sp(sp,13)} = [npoints{CONNECT_sp(sp,13)} [sp;1]];
 CONNECT_sp(sp,14)= CONNECT_sp(sp,13) + 1;                        
npoints{CONNECT_sp(sp,14)} = [npoints{CONNECT_sp(sp,14)} [sp;1]];
 CONNECT_sp(sp,15)= CONNECT_sp(sp,14) + 1;                        
npoints{CONNECT_sp(sp,15)} = [npoints{CONNECT_sp(sp,15)} [sp;1]];
 CONNECT_sp(sp,16)= CONNECT_sp(sp,15) + 1;                        
npoints{CONNECT_sp(sp,16)} = [npoints{CONNECT_sp(sp,16)} [sp;1]];

 % Node index of data position including particle and corners
%  CONNECT_MLS{sp,1}(1:4) = CONNECT_sp(sp,1:4);
 CONNECT_MLS{sp,1}(1:16) = CONNECT_sp(sp,1:16);

%  % Particle-cell interaction
%  CONNECT_sp_G(sp,1) = spElems(sp)-(NN(1)-1)-1;     
%  cpoints{CONNECT_sp_G(sp,1)} = [cpoints{CONNECT_sp_G(sp,1)} [sp;1]];
%  CONNECT_sp_G(sp,2) = spElems(sp)-(NN(1)-1);       
%  cpoints{CONNECT_sp_G(sp,2)} = [cpoints{CONNECT_sp_G(sp,2)} [sp;1]];
%  CONNECT_sp_G(sp,3) = spElems(sp)-(NN(1)-1)+1;     
%  cpoints{CONNECT_sp_G(sp,3)} = [cpoints{CONNECT_sp_G(sp,3)} [sp;1]];
%  CONNECT_sp_G(sp,4) = spElems(sp)-1;               
%  cpoints{CONNECT_sp_G(sp,4)} = [cpoints{CONNECT_sp_G(sp,4)} [sp;1]];
%  CONNECT_sp_G(sp,5) = spElems(sp);                 
%  cpoints{CONNECT_sp_G(sp,5)} = [cpoints{CONNECT_sp_G(sp,5)} [sp;1]];
%  CONNECT_sp_G(sp,6) = spElems(sp)+1;               
%  cpoints{CONNECT_sp_G(sp,6)} = [cpoints{CONNECT_sp_G(sp,6)} [sp;1]];
%  CONNECT_sp_G(sp,7) = spElems(sp)+(NN(1)-1)-1;     
%  cpoints{CONNECT_sp_G(sp,7)} = [cpoints{CONNECT_sp_G(sp,7)} [sp;1]];
%  CONNECT_sp_G(sp,8) = spElems(sp)+(NN(1)-1);       
%  cpoints{CONNECT_sp_G(sp,8)} = [cpoints{CONNECT_sp_G(sp,8)} [sp;1]];
%  CONNECT_sp_G(sp,9) = spElems(sp)+(NN(1)-1)+1;     
%  cpoints{CONNECT_sp_G(sp,9)} = [cpoints{CONNECT_sp_G(sp,9)} [sp;1]];
% 
%   % Cell index of data position including particle and corners
%  CONNECT_MLS_G{sp,1}(1:9) = CONNECT_sp_G(sp,1:9);
 end
 
 %% Corner interaction
 for sp = 1:spCount
 % Corner 1-node interaction 
 [CONNECT_MLS,npoints]=FindNodeInteractingParticle...
    (sp,spElems_corner(sp,1),npoints,NN,2,CONNECT_MLS);
 % Corner 2-node interaction
 [CONNECT_MLS,npoints]=FindNodeInteractingParticle...
    (sp,spElems_corner(sp,2),npoints,NN,3,CONNECT_MLS);
 % Corner 3-node interaction
 [CONNECT_MLS,npoints]=FindNodeInteractingParticle...
    (sp,spElems_corner(sp,3),npoints,NN,4,CONNECT_MLS);
 % Corner 4-node interaction
 [CONNECT_MLS,npoints]=FindNodeInteractingParticle...
    (sp,spElems_corner(sp,4),npoints,NN,5,CONNECT_MLS);
%  
%  % Corner 1-cell interaction
%  [CONNECT_MLS_G,cpoints]=FindCellInteractingParticle...
%     (sp,spElems_corner(sp,1),cpoints,NN,2,CONNECT_MLS_G);
%  % Corner 2-cell interaction
% [CONNECT_MLS_G,cpoints]=FindCellInteractingParticle...
%     (sp,spElems_corner(sp,2),cpoints,NN,3,CONNECT_MLS_G);
%  % Corner 3-cell interaction
% [CONNECT_MLS_G,cpoints]=FindCellInteractingParticle...
%     (sp,spElems_corner(sp,3),cpoints,NN,4,CONNECT_MLS_G);
%  % Corner 4-cell interaction
% [CONNECT_MLS_G,cpoints]=FindCellInteractingParticle...
%     (sp,spElems_corner(sp,4),cpoints,NN,5,CONNECT_MLS_G);
 end
 
%  
%   % Node index for centroid
%  for c = 1:cellCount
%  cElems(c) = ceil(LOCC(c,1)/le(1))+(NN(1)-1)*(fix(LOCC(c,2)/le(2)));
%  CONNECT_TEMP_C(c,1) = cElems(c) + ceil(cElems(c)/(NN(1)-1)) - 1;
%  CONNECT_TEMP_C(c,2) = CONNECT_TEMP_C(c,1)+1; 
%  CONNECT_TEMP_C(c,3) = CONNECT_TEMP_C(c,2)+NN(1); 
%  CONNECT_TEMP_C(c,4) = CONNECT_TEMP_C(c,1)+NN(1);
%  CONNECT_C{c}        = [CONNECT_TEMP_C(c,1) CONNECT_TEMP_C(c,2) CONNECT_TEMP_C(c,3) CONNECT_TEMP_C(c,4)];
%  end
 
 % List of particles index inside each cell "c"
 for c =1:cellCount
     id_sp = find(spElems==c);
     mspoints{c}=id_sp;
 end
 
%   %% Linear for mapping from centroids to node
%  for c = 1:cellCount
%      for i=1:4
%      npid = CONNECT_TEMP_C(c,i);
%      % Compute the shape functions and gradient of the shape functions
%     [Nc_local(c,i),dNcx_local(c,i),dNcy_local(c,i)]=linearshape(LOCC(c,:),LOC(npid,:),le(1),le(2));
%     N_c{c}(i) = Nc_local(c,i);    
%     dN_c{c}(1,i) = dNcx_local(c,i);
%     dN_c{c}(2,i) = dNcy_local(c,i);
%      end
%  end 
 
 %% ---------------------------------------------------------End CPDI interpolator
%  %% Determine the cell boundary of object
%     Ncorner_boundary = length(boundary_corner(1,:));
%         
% %  spElems_boundary = zeros(1,5*Ncorner_boundary);
%   spElems_boundary = zeros(1,Ncorner_boundary);
% 
%  for i = 1:Ncorner_boundary
%      index_corner = boundary_corner(1:2,i);
%      particle_index = index_corner(1);
%      corner_index = index_corner(2);
%      spElems_boundary(i) = ceil(x_data{particle_index,corner_index}(1)/le(1))+(NN(1)-1)*(fix(x_data{particle_index,corner_index}(2)/le(2)));
% 
% %      for j = 1:5
% %      spElems_boundary(5*(i-1)+j) = ceil(x_data{particle_index,j}(1)/le(1))+(NN(1)-1)*(fix(x_data{particle_index,j}(2)/le(2)));
% %      end
%      
%  end
%  spElems_boundary = unique(spElems_boundary);
% 
%  %% Compute cell boundary node
%  CellBoundary_node = [];
%  NCellBoundary = length(spElems_boundary);
%  
%  for i = 1:NCellBoundary
%      cpid = spElems_boundary(i);
%      
%      CellBoundary_node = [CellBoundary_node CONNECT_C{cpid}];
%  end
% %  CellBoundary_node = [CellBoundary_node node_boundary];
%  CellBoundary_node = unique(CellBoundary_node);
 
 
 %% Build IMLS_CPDI interpolator
 Elems_CPDI  = zeros(spCount,5);            % Elements store the particles and corners
  
 for sp = 1:spCount
    Elems_CPDI(sp,:) = [spElems(sp) spElems_corner(sp,1) spElems_corner(sp,2) spElems_corner(sp,3) spElems_corner(sp,4)];
    
%      % Cell index of data position including particle and corners
%      for i = 1:9
%      CONNECT_MLS_G{sp,2} = [CONNECT_MLS_G{sp,2} CONNECT_TEMP_G(sp,i)];
%      CONNECT_MLS_G{sp,3} = [CONNECT_MLS_G{sp,3} CONNECT_TEMP_G(sp,i+9)];
%      CONNECT_MLS_G{sp,4} = [CONNECT_MLS_G{sp,4} CONNECT_TEMP_G(sp,i+18)];
%      CONNECT_MLS_G{sp,5} = [CONNECT_MLS_G{sp,5} CONNECT_TEMP_G(sp,i+27)];
%      end
 end

%  %% Quadratic-Bspline for mapping from particle to centroids using IMLS
%  [N_g] = CellMLS(spElems_boundary,spCount,cellCount,cpoints,x_data,LOCC,le,Elems_CPDI,CONNECT_MLS_G);

 %% Cubic spline for mapping from particle to nodes using IMLS 
 [k,node_boundary,color_node,N_MLS] = NodeMLS(spCount,nodeCount,npoints,x_data,LOC,le,Elems_CPDI,CONNECT_MLS);

end