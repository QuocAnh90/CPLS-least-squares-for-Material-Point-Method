function [N,dN,N_MLS,N_g,N_c,dN_c,CONNECT,CONNECT_G,CONNECT_C,CONNECT_N...
,spElems,mspoints,NODES,GAUSS,NODES_G,NODES_N] = Compute_Interpolator_IMLS_MPM(spCount,cellCount,nodeCount,x_sp,le,NN,LOC,LOCC)

%% Output
% Mapping from particles to nodes using linear basis N{p}(i); dN{p}(1:2,i)
% Mapping from particles to nodes using moving least square N_MLS{p}(i)
% Mapping from particles to centroids using quadratic Bspline N_g{p}(c)
% Mapping from centroids to nodes using linear basis Nc{p}(i); dNc{p}(1:2,i)


%% Parameter generations
spElems         = zeros(spCount,1);              % index of elements where stores particles
CONNECT_TEMP    = zeros(spCount,4);              % node 1=leftdown 2=righdown 3=rightup 4= leftup
CONNECT_TEMP_G  = zeros(spCount,9);              % gauss number from left to right from bottom to top
CONNECT_TEMP_C  = zeros(cellCount,4);
CONNECT_TEMP_N  = zeros(spCount,16);  
NODES           = 4 * ones(spCount,1);           % Number of interation nodes for each particle
NODES_N         = 16 * ones(spCount,1);
GAUSS           = 9 * ones(spCount,1);
cElems          = zeros(cellCount,1);
NODES_G         = 4 * ones(cellCount,1);          
cpoints         = cell(cellCount,1);
cpoints1        = cell(nodeCount,1);
npoints         = cell(nodeCount,1);
npoints1        = cell(nodeCount,1);
CONNECT_N       = cell(spCount,1);

% Funtion for mapping from particles to nodes using linear basis
N_local         = zeros(spCount,4);              % Value of shape function
dNx_local       = zeros(spCount,4);              % Value of gradient of shape function
dNy_local       = zeros(spCount,4); 

% Funtion for mapping from particles to nodes using cubic spline
Nn_local        = zeros(spCount,16);              % Value of shape function

% Funtion for mapping from particles to nodes using moving least square
% Least square matrix
A               = cell(nodeCount,1);
dAx             = cell(nodeCount,1);
dAy             = cell(nodeCount,1);

% Function for mapping from particles to centroids using quadratic Bspline
Ng_local         = zeros(spCount,9);              % Value of shape function
N_g              = cell(spCount,1);
% Least square matrix
Ac               = cell(cellCount,1);

% Funtion for mapping from centroids to nodes using linear basis
Nc_local         = zeros(cellCount,4);              % Value of shape function
dNcx_local       = zeros(cellCount,4);              % Value of gradient of shape function
dNcy_local       = zeros(cellCount,4); 

% Order of MLS
r = 3;
r_n = r * ones(nodeCount,1);
r_c = r * ones(nodeCount,1);

for i = 1:nodeCount
    A{i} = zeros(r_n(i),r_n(i));
    dAx{i} = zeros(r_n(i),r_n(i));
    dAy{i} = zeros(r_n(i),r_n(i));
end

for c = 1:cellCount
    Ac{c} = zeros(r_c(c),r_c(c));
end
 
%% Compute the index of the node for the particles 
% and the index of the centroids to particles

 for p = 1:spCount
%  spElems(p) = ceil(x_sp(p,1)/le(1))+(NN(1)-1)*(fix(x_sp(p,2)/le(2)));   % compute vector store index elements 
 spElems(p) = floor(x_sp(p,1)/le(1)+1)+(NN(1)-1)*(floor(x_sp(p,2)/le(2)));
 
 % Node index for particle using linear basis
 CONNECT_TEMP(p,1) = spElems(p) + ceil(spElems(p)/(NN(1)-1)) - 1;       npoints1{CONNECT_TEMP(p,1)} = [npoints1{CONNECT_TEMP(p,1)} p];
 CONNECT_TEMP(p,2) = CONNECT_TEMP(p,1)+1;                               npoints1{CONNECT_TEMP(p,2)} = [npoints1{CONNECT_TEMP(p,2)} p];
 CONNECT_TEMP(p,3) = CONNECT_TEMP(p,2)+NN(1);                           npoints1{CONNECT_TEMP(p,3)} = [npoints1{CONNECT_TEMP(p,3)} p];
 CONNECT_TEMP(p,4) = CONNECT_TEMP(p,1)+NN(1);                           npoints1{CONNECT_TEMP(p,4)} = [npoints1{CONNECT_TEMP(p,4)} p];
 CONNECT{p}        = [CONNECT_TEMP(p,1) CONNECT_TEMP(p,2) CONNECT_TEMP(p,3) CONNECT_TEMP(p,4)];
 
 % Node index for particle using cubic spline
 CONNECT_TEMP_N(p,6) = spElems(p) + ceil(spElems(p)/(NN(1)-1)) - 1;     npoints{CONNECT_TEMP_N(p,6)} = [npoints{CONNECT_TEMP_N(p,6)} p];
 CONNECT_TEMP_N(p,7) = CONNECT_TEMP_N(p,6) + 1;                         npoints{CONNECT_TEMP_N(p,7)} = [npoints{CONNECT_TEMP_N(p,7)} p];
 CONNECT_TEMP_N(p,11)= CONNECT_TEMP_N(p,7) + NN(1);                     npoints{CONNECT_TEMP_N(p,11)} = [npoints{CONNECT_TEMP_N(p,11)} p];
 CONNECT_TEMP_N(p,10)= CONNECT_TEMP_N(p,6) + NN(1);                     npoints{CONNECT_TEMP_N(p,10)} = [npoints{CONNECT_TEMP_N(p,10)} p];
 CONNECT_TEMP_N(p,1) = CONNECT_TEMP_N(p,6) - NN(1) - 1;                 npoints{CONNECT_TEMP_N(p,1)} = [npoints{CONNECT_TEMP_N(p,1)} p];
 CONNECT_TEMP_N(p,2) = CONNECT_TEMP_N(p,1) + 1;                         npoints{CONNECT_TEMP_N(p,2)} = [npoints{CONNECT_TEMP_N(p,2)} p];
 CONNECT_TEMP_N(p,3) = CONNECT_TEMP_N(p,2) + 1;                         npoints{CONNECT_TEMP_N(p,3)} = [npoints{CONNECT_TEMP_N(p,3)} p];
 CONNECT_TEMP_N(p,4) = CONNECT_TEMP_N(p,3) + 1;                         npoints{CONNECT_TEMP_N(p,4)} = [npoints{CONNECT_TEMP_N(p,4)} p];
 CONNECT_TEMP_N(p,5) = CONNECT_TEMP_N(p,6) - 1;                         npoints{CONNECT_TEMP_N(p,5)} = [npoints{CONNECT_TEMP_N(p,5)} p];
 CONNECT_TEMP_N(p,8) = CONNECT_TEMP_N(p,7) + 1;                         npoints{CONNECT_TEMP_N(p,8)} = [npoints{CONNECT_TEMP_N(p,8)} p];
 CONNECT_TEMP_N(p,9) = CONNECT_TEMP_N(p,5) + NN(1);                     npoints{CONNECT_TEMP_N(p,9)} = [npoints{CONNECT_TEMP_N(p,9)} p];
 CONNECT_TEMP_N(p,12)= CONNECT_TEMP_N(p,11) + 1;                        npoints{CONNECT_TEMP_N(p,12)} = [npoints{CONNECT_TEMP_N(p,12)} p];
 CONNECT_TEMP_N(p,13)= CONNECT_TEMP_N(p,9) + NN(1);                     npoints{CONNECT_TEMP_N(p,13)} = [npoints{CONNECT_TEMP_N(p,13)} p];
 CONNECT_TEMP_N(p,14)= CONNECT_TEMP_N(p,13) + 1;                        npoints{CONNECT_TEMP_N(p,14)} = [npoints{CONNECT_TEMP_N(p,14)} p];
 CONNECT_TEMP_N(p,15)= CONNECT_TEMP_N(p,14) + 1;                        npoints{CONNECT_TEMP_N(p,15)} = [npoints{CONNECT_TEMP_N(p,15)} p];
 CONNECT_TEMP_N(p,16)= CONNECT_TEMP_N(p,15) + 1;                        npoints{CONNECT_TEMP_N(p,16)} = [npoints{CONNECT_TEMP_N(p,16)} p];

 for i=1:16
 CONNECT_N{p} = [CONNECT_N{p} CONNECT_TEMP_N(p,i)];
 end
 
 % Centroid index for particle
 CONNECT_TEMP_G(p,1) =  spElems(p)-(NN(1)-1)-1;     cpoints{CONNECT_TEMP_G(p,1)} = [cpoints{CONNECT_TEMP_G(p,1)} p];
 CONNECT_TEMP_G(p,2) =  spElems(p)-(NN(1)-1);       cpoints{CONNECT_TEMP_G(p,2)} = [cpoints{CONNECT_TEMP_G(p,2)} p];
 CONNECT_TEMP_G(p,3) =  spElems(p)-(NN(1)-1)+1;     cpoints{CONNECT_TEMP_G(p,3)} = [cpoints{CONNECT_TEMP_G(p,3)} p];
 CONNECT_TEMP_G(p,4) =  spElems(p)-1;               cpoints{CONNECT_TEMP_G(p,4)} = [cpoints{CONNECT_TEMP_G(p,4)} p];
 CONNECT_TEMP_G(p,5) =  spElems(p);                 cpoints{CONNECT_TEMP_G(p,5)} = [cpoints{CONNECT_TEMP_G(p,5)} p];
 CONNECT_TEMP_G(p,6) =  spElems(p)+1;               cpoints{CONNECT_TEMP_G(p,6)} = [cpoints{CONNECT_TEMP_G(p,6)} p];
 CONNECT_TEMP_G(p,7) =  spElems(p)+(NN(1)-1)-1;     cpoints{CONNECT_TEMP_G(p,7)} = [cpoints{CONNECT_TEMP_G(p,7)} p];
 CONNECT_TEMP_G(p,8) =  spElems(p)+(NN(1)-1);       cpoints{CONNECT_TEMP_G(p,8)} = [cpoints{CONNECT_TEMP_G(p,8)} p];
 CONNECT_TEMP_G(p,9) =  spElems(p)+(NN(1)-1)+1;     cpoints{CONNECT_TEMP_G(p,9)} = [cpoints{CONNECT_TEMP_G(p,9)} p];
 CONNECT_G{p} = [CONNECT_TEMP_G(p,1) CONNECT_TEMP_G(p,2) CONNECT_TEMP_G(p,3) CONNECT_TEMP_G(p,4) CONNECT_TEMP_G(p,5) CONNECT_TEMP_G(p,6) CONNECT_TEMP_G(p,7) CONNECT_TEMP_G(p,8) CONNECT_TEMP_G(p,9)];
 end
 
 % Node index for centroid
 for c = 1:cellCount
 cElems(c) = ceil(LOCC(c,1)/le(1))+(NN(1)-1)*(fix(LOCC(c,2)/le(2)));
 CONNECT_TEMP_C(c,1) = cElems(c) + ceil(cElems(c)/(NN(1)-1)) - 1;
 CONNECT_TEMP_C(c,2) = CONNECT_TEMP_C(c,1)+1; 
 CONNECT_TEMP_C(c,3) = CONNECT_TEMP_C(c,2)+NN(1); 
 CONNECT_TEMP_C(c,4) = CONNECT_TEMP_C(c,1)+NN(1);
 CONNECT_C{c}        = [CONNECT_TEMP_C(c,1) CONNECT_TEMP_C(c,2) CONNECT_TEMP_C(c,3) CONNECT_TEMP_C(c,4)];
 end
     
%% Compute mspoints: index of particles in each element (active element)
 for c =1:cellCount
     id_p = find(spElems==c);
     mspoints{c}=id_p;
 end
 
%% Linear basis function from particle to nodes
 for p = 1:spCount
for i = 1:NODES(p)
    npid = CONNECT_TEMP(p,i);
     % Compute the shape functions and gradient of the shape functions
    [N_local(p,i),dNx_local(p,i),dNy_local(p,i)]=linearshape(x_sp(p,1:2),LOC(npid,:),le(1),le(2));
    N{p}(i) = N_local(p,i);    
    dN{p}(1,i) = dNx_local(p,i);
    dN{p}(2,i) = dNy_local(p,i);
end
end

 
% %% Quadratic-Bspline for mapping from particle to centroids
%  for p = 1:spCount
% for c = 1:GAUSS(p)
%     cpid = CONNECT_TEMP_G(p,c);
%      % Compute the shape functions and gradient of the shape functions
%     [Ng_local(p,c)]=Quadratic_Bspline(x_sp(p,1:2),LOCC(cpid,:),le(1),le(2));
% end
%  end
%  
% % Least square function from particle to centroids
% % Compute A
% for p = 1:spCount
%      if r==1
%          Pxy = 1;
%      elseif r==3
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];    
%      elseif r==6
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
%      end
%          
%      PPxy = Pxy * Pxy';
%      
%      for c = 1:GAUSS(p)
%          cpid = CONNECT_TEMP_G(p,c);
%          Ac{cpid} = Ac{cpid} + Ng_local(p,c) * PPxy;
%      end
% end
% 
% % % Check rank D
% % for c=1:cellCount
% % %     if rank(A{c})<6
% %     if rcond(A{c})<=9.999999e-10
% %         r_c(c) = 3;
% %     end
% % end
%  
% % Compute shape function and gradients
%  for p = 1:spCount
%      if r==1
%          Pxy = 1;
%      elseif r==3
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];   
%      elseif r==6
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
%      end
%      
%      for c = 1:GAUSS(p)
%          cpid = CONNECT_TEMP_G(p,c);
%                 
%          if r==1
%          pc = [1];
%          elseif r==3
%          pc = [1 ; LOCC(cpid,1) ; LOCC(cpid,2)];   
%          elseif r==6
%          pc = [1 ; LOCC(cpid,1) ; LOCC(cpid,2) ; LOCC(cpid,1)^2 ; LOCC(cpid,1)*LOCC(cpid,2) ; LOCC(cpid,2)^2];   
%          end
%          
%          if r_c(cpid)==6
%          
%          rc = Ac{cpid} \ pc;
%          N_g{p}(c) = rc' * Ng_local(p,c) * Pxy;  
%          
%          elseif r_c(cpid)==3
%          rc = Ac{cpid}(1:3,1:3) \ pc(1:3);
%          N_g{p}(c) = rc' * Ng_local(p,c) * Pxy(1:3);  
%          end
%          
%      end
%  end
 
 %% Quadratic-Bspline for mapping from particle to centroids using IMLS
 for c = 1:cellCount
%      if isempty(mspoints{c})==1
%          continue
%      end
     nv = x_sp(cpoints{c},:);
     nv = nv';
     Np = length(cpoints{c});
     won = ones(1,Np);
     w = zeros(1,Np);
     for i = 1:Np
         pid = cpoints{c}(i);
         Ng_local(pid,c) = Quadratic_Bspline(x_sp(pid,1:2),LOCC(c,:),le(1),le(2));
         w(i) = Ng_local(pid,c);
     end
     p = [won;nv];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=zeros(3,Np);
q(1,1:Np)=1;
gama=zeros(1,3);
bita=zeros(1,3);
Rfa=zeros(1,3);

p1xi=p(1,1:Np);
q1xi=q(1,1:Np);

p2xi=p(2,1:Np);
gama1=0;
for i=1:Np
    gama1=gama1+w(i)*q1xi(i)^2;
end
gama(1)=gama1;
bita21=0;
for i=1:Np
    bita21=bita21+w(i)*p2xi(i)*q1xi(i);
end
bita(1)=bita21;
if gama1==0
    Rfa21=0;
else
    Rfa21=bita21/gama1;
end

Rfa(1)=Rfa21;


for i=1:Np
    q(2,i)=p2xi(i)-Rfa21*q1xi(i);
end
q2xi=q(2,1:Np);

p3xi=p(3,1:Np);
gama2=0;
bita31=0;
bita32=0;
for i=1:Np
    gama2=gama2+w(i)*q2xi(i)^2;
    bita31=bita31+w(i)*p3xi(i)*q1xi(i);
    bita32=bita32+w(i)*p3xi(i)*q2xi(i);
end
gama(2)=gama2;
bita(2)=bita31;
bita(3)=bita32;

if gama1==0
    Rfa31=0;
else
    Rfa31=bita31/gama1;
end

if gama2==0
    Rfa32=0;
else
    Rfa32=bita32/gama2;
end

Rfa(2)=Rfa31;
Rfa(3)=Rfa32;
for i=1:Np
    q(3,i)=p3xi(i)-Rfa31*q1xi(i)-Rfa32*q2xi(i);
end

q3xi=q(3,1:Np);

gama3=0;
for i=1:Np
    gama3=gama3+w(i)*q3xi(i)^2;
end
gama(3)=gama3;

gposx=LOCC(c,1);
gposy=LOCC(c,2);
qq=zeros(1,3);
qq(1)=1;
qq(2)=gposx-Rfa(1);
qq(3)=gposy-(Rfa(2)*qq(1)+Rfa(3)*qq(2));
cc=zeros(3,Np);
for j=1:3
    for i=1:Np
        if gama(j)==0
            cc(j,i)=0;
        else
        cc(j,i)=(qq(j)*q(j,i))/gama(j);
        end
    end
end
phi=zeros(1,Np);
for i=1:Np
    cji=0;
    for j=1:3
        cji=cji+cc(j,i);
    end
    phi(i)=w(i)*cji;
end

     for i = 1:Np
         pid = cpoints{c}(i);
     id_c = find(CONNECT_TEMP_G(pid,:)==c);
     N_g{pid}(id_c) = phi(i);
     end
 end

 %% Linear for mapping from centroids to node
 for c = 1:cellCount
     for i=1:NODES_G(c)
     npid = CONNECT_TEMP_C(c,i);
     % Compute the shape functions and gradient of the shape functions
    [Nc_local(c,i),dNcx_local(c,i),dNcy_local(c,i)]=linearshape(LOCC(c,:),LOC(npid,:),le(1),le(2));
    N_c{c}(i) = Nc_local(c,i);    
    dN_c{c}(1,i) = dNcx_local(c,i);
    dN_c{c}(2,i) = dNcy_local(c,i);
     end
 end 

% %% Least square function from particle to nodes using linear shape
%  for p = 1:spCount
% for i = 1:NODES(p)
%     npid = CONNECT_TEMP(p,i);
%      % Compute the shape functions and gradient of the shape functions
%     [N_local(p,i),~,~]=linearshape(x_sp(p,1:2),LOC(npid,:),le(1),le(2));
% end
%  end
% 
% % Compute A
% for p = 1:spCount
%      if r==1
%          Pxy = 1;
%      elseif r==3
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];    
%      elseif r==6
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
%      end
%          
%      PPxy = Pxy * Pxy';
%      
%      for i = 1:NODES(p)
%          npid = CONNECT{p}(i);
%          A{npid} = A{npid} + N_local(p,i) * PPxy;  
%      end
% end
% 
% % % Check rank D
% % for n=1:nodeCount
% % %     if rank(A{n})<6
% %     if rcond(A{n})<=9.999999e-10
% %         r_n(n) = 3;
% %     end
% % end
% 
% % Compute shape function and gradients
%  for p = 1:spCount
%      if r==1
%          Pxy = 1;
%      elseif r==3
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];   
%      elseif r==6
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
%      end
%      
%      for i = 1:NODES(p)
%          npid = CONNECT{p}(i);
%          if r==1
%          pn = 1;
%          elseif r==3
%          pn = [1 ; LOC(npid,1) ; LOC(npid,2)];   
%          elseif r==6
%          pn = [1 ; LOC(npid,1) ; LOC(npid,2) ; LOC(npid,1)^2 ; LOC(npid,1)*LOC(npid,1) ; LOC(npid,1)^2];   
%          end
%          
%          if r_n(npid)==6
%          rn = A{npid} \ pn;
%          N_MLS{p}(i) = rn' * N_local(p,i) * Pxy;
%          
%          elseif r_n(npid)==3
%          rn = A{npid}(1:3,1:3) \ pn(1:3);
%          N_MLS{p}(i) = rn' * N_local(p,i) * Pxy(1:3);       
%          end        
%      end
%  end

% %% Least square function from particle to nodes using cubic spline
%  for p = 1:spCount
% for i = 1:NODES_N(p)
%     npid = CONNECT_TEMP_N(p,i);
%      % Compute the shape functions and gradient of the shape functions
%     [Nn_local(p,i)]=Cubic_Bspline(x_sp(p,1:2),LOC(npid,:),le(1),le(2));
% end
%  end
% 
% % Compute A
% for p = 1:spCount
%      if r==1
%          Pxy = 1;
%      elseif r==3
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];    
%      elseif r==6
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
%      end
%          
%      PPxy = Pxy * Pxy';
%      
%      for i = 1:NODES_N(p)
%          npid = CONNECT_N{p}(i);
%          A{npid} = A{npid} + Nn_local(p,i) * PPxy;  
%      end
% end
% 
% % % Check rank D
% % for n=1:nodeCount
% % %     if rank(A{n})<6
% %     if rcond(A{n})<=9.999999e-10
% %         r_n(n) = 3;
% %     end
% % end
% 
% % Compute shape function and gradients
%  for p = 1:spCount
%      if r==1
%          Pxy = 1;
%      elseif r==3
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];   
%      elseif r==6
%          Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
%      end
%      
%      for i = 1:NODES_N(p)
%          npid = CONNECT_N{p}(i);
%          if r==1
%          pn = 1;
%          elseif r==3
%          pn = [1 ; LOC(npid,1) ; LOC(npid,2)];   
%          elseif r==6
%          pn = [1 ; LOC(npid,1) ; LOC(npid,2) ; LOC(npid,1)^2 ; LOC(npid,1)*LOC(npid,1) ; LOC(npid,1)^2];   
%          end
%          
%          if r_n(npid)==6
%          rn = A{npid} \ pn;
%          N_MLS{p}(i) = rn' * Nn_local(p,i) * Pxy;
%          
%          elseif r_n(npid)==3
%          rn = A{npid}(1:3,1:3) \ pn(1:3);
%          N_MLS{p}(i) = rn' * Nn_local(p,i) * Pxy(1:3);
%            
%          end
%          
%      end
%  end

 %% Linear for mapping from particles to nodes using IMLS with linear basis
 for n = 1:nodeCount
     if isempty(npoints1{n})==1
         continue
     end
     cv = x_sp(npoints1{n},:);
     cv = cv';
     Np = length(npoints1{n});
     won = ones(1,Np);
     wn = zeros(1,Np);
     for i = 1:Np
         pid = npoints1{n}(i);
         N_local(pid,n) = linearshape(x_sp(pid,1:2),LOC(n,:),le(1),le(2));
%          Nn_local(pid,n) = Cubic_Bspline(x_sp(pid,1:2),LOC(n,:),le(1),le(2));
         wn(i) = N_local(pid,n);
     end
     p = [won;cv];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=zeros(3,Np);
q(1,1:Np)=1;
gama=zeros(1,3);
bita=zeros(1,3);
Rfa=zeros(1,3);

p1xi=p(1,1:Np);
q1xi=q(1,1:Np);

p2xi=p(2,1:Np);
gama1=0;
for i=1:Np
    gama1=gama1+wn(i)*q1xi(i)^2;
end
gama(1)=gama1;
bita21=0;
for i=1:Np
    bita21=bita21+wn(i)*p2xi(i)*q1xi(i);
end
bita(1)=bita21;

if gama1==0
    Rfa21=0;
else
    Rfa21=bita21/gama1;
end

Rfa(1)=Rfa21;

for i=1:Np
    q(2,i)=p2xi(i)-Rfa21*q1xi(i);
end
q2xi=q(2,1:Np);

p3xi=p(3,1:Np);
gama2=0;
bita31=0;
bita32=0;
for i=1:Np
    gama2=gama2+wn(i)*q2xi(i)^2;
    bita31=bita31+wn(i)*p3xi(i)*q1xi(i);
    bita32=bita32+wn(i)*p3xi(i)*q2xi(i);
end
gama(2)=gama2;
bita(2)=bita31;
bita(3)=bita32;

if gama1==0
    Rfa31=0;
else
    Rfa31=bita31/gama1;
end

if gama2==0
    Rfa32=0;
else
    Rfa32=bita32/gama2;
end

Rfa(2)=Rfa31;
Rfa(3)=Rfa32;
for i=1:Np
    q(3,i)=p3xi(i)-Rfa31*q1xi(i)-Rfa32*q2xi(i);
end

q3xi=q(3,1:Np);

gama3=0;
for i=1:Np
    gama3=gama3+wn(i)*q3xi(i)^2;
end
gama(3)=gama3;

gposx=LOC(n,1);
gposy=LOC(n,2);

qq=zeros(1,3);
qq(1)=1;
qq(2)=gposx-Rfa(1);
qq(3)=gposy-(Rfa(2)*qq(1)+Rfa(3)*qq(2));
cc=zeros(3,Np);
for j=1:3
    for i=1:Np
        if gama(j)==0
            cc(j,i)=0;
        else
            cc(j,i)=(qq(j)*q(j,i))/gama(j);
        end
    end
end
phi=zeros(1,Np);
for i=1:Np
    cji=0;
    for j=1:3
        cji=cji+cc(j,i);
    end
    phi(i)=wn(i)*cji;
end

     for i = 1:Np
         pid = npoints1{n}(i);
     id_n = find(CONNECT_TEMP(pid,:)==n);
%      id_n = find(CONNECT_TEMP_N(pid,:)==n);
     N_MLS{pid}(id_n) = phi(i);
     end
 end

% %% Cubic spline for mapping from particle to nodes using IMLS
%  for n = 1:nodeCount
%      if isempty(npoints{n})==1
%          continue
%      end
%      cv = x_sp(npoints{n},:);
%      cv = cv';
%      Np = length(npoints{n});
%      won = ones(1,Np);
%      wn = zeros(1,Np);
%      for i = 1:Np
%          pid = npoints{n}(i);
% %          N_local(pid,n) = linearshape(x_sp(pid,1:2),LOC(n,:),le(1),le(2));
%          Nn_local(pid,n) = Cubic_Bspline(x_sp(pid,1:2),LOC(n,:),le(1),le(2));
%          wn(i) = Nn_local(pid,n);
%      end
%      p = [won;cv];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q=zeros(3,Np);
% q(1,1:Np)=1;
% gama=zeros(1,3);
% bita=zeros(1,3);
% Rfa=zeros(1,3);
% 
% p1xi=p(1,1:Np);
% q1xi=q(1,1:Np);
% 
% p2xi=p(2,1:Np);
% gama1=0;
% for i=1:Np
%     gama1=gama1+wn(i)*q1xi(i)^2;
% end
% gama(1)=gama1;
% bita21=0;
% for i=1:Np
%     bita21=bita21+wn(i)*p2xi(i)*q1xi(i);
% end
% bita(1)=bita21;
% if gama1==0
%     Rfa21=0;
% else
%     Rfa21=bita21/gama1;
% end
% Rfa(1)=Rfa21;
% 
% for i=1:Np
%     q(2,i)=p2xi(i)-Rfa21*q1xi(i);
% end
% q2xi=q(2,1:Np);
% 
% p3xi=p(3,1:Np);
% gama2=0;
% bita31=0;
% bita32=0;
% for i=1:Np
%     gama2=gama2+wn(i)*q2xi(i)^2;
%     bita31=bita31+wn(i)*p3xi(i)*q1xi(i);
%     bita32=bita32+wn(i)*p3xi(i)*q2xi(i);
% end
% gama(2)=gama2;
% bita(2)=bita31;
% bita(3)=bita32;
% if gama1==0
%     Rfa31=0;
% else
%     Rfa31=bita31/gama1;
% end
% 
% if gama2==0
%     Rfa32=0;
% else
%     Rfa32=bita32/gama2;
% end
% 
% Rfa(2)=Rfa31;
% Rfa(3)=Rfa32;
% for i=1:Np
%     q(3,i)=p3xi(i)-Rfa31*q1xi(i)-Rfa32*q2xi(i);
% end
% 
% q3xi=q(3,1:Np);
% 
% gama3=0;
% for i=1:Np
%     gama3=gama3+wn(i)*q3xi(i)^2;
% end
% gama(3)=gama3;
% 
% gposx=LOC(n,1);
% gposy=LOC(n,2);
% qq=zeros(1,3);
% qq(1)=1;
% qq(2)=gposx-Rfa(1);
% qq(3)=gposy-(Rfa(2)*qq(1)+Rfa(3)*qq(2));
% cc=zeros(3,Np);
% for j=1:3
%     for i=1:Np
%         if gama(j)==0
%         cc(j,i)=0;
%         else
%         cc(j,i)=(qq(j)*q(j,i))/gama(j);
%         end
%     end
% end
% phi=zeros(1,Np);
% for i=1:Np
%     cji=0;
%     for j=1:3
%         cji=cji+cc(j,i);
%     end
%     phi(i)=wn(i)*cji;
% end
% 
%      for i = 1:Np
%          pid = npoints{n}(i);
% %      id_n = find(CONNECT_TEMP(pid,:)==n);
%      id_n = find(CONNECT_TEMP_N(pid,:)==n);
%      N_MLS{pid}(id_n) = phi(i);
%      end
%  end
 
 test = 1;
