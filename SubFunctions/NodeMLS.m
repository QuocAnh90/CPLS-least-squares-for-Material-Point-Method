function [k,node_boundary,color_node,N_MLS] = NodeMLS...
    (spCount,nodeCount,npoints,x_data,LOC,le,Elems_CPDI,CONNECT_MLS)

%% Cubic spline for mapping from particle to nodes using IMLS
  node_boundary  = [];
  color_node     = zeros(nodeCount,1);
  k              = zeros(nodeCount,1);
  N_MLS          = cell(spCount,5);
%   Number_cell    = zeros(nodeCount,1);

  
 for n = 1:nodeCount
     if isempty(npoints{n})==1
         continue
     end

     Np = length(npoints{n}(1,:));  
%     % Compute MLS shape function using node boundaries
%     Ncell = [];          % Node contain cells with particles
%     
%     for i = 1:Np
%      pid = npoints{n}(1,i);
%      corner = npoints{n}(2,i);
%      id_n = find(CONNECT_MLS{pid,corner}==n);
%      
%      if id_n == 6 || id_n == 7 || id_n == 10 || id_n == 11 
%      Ncell = [Ncell Elems_CPDI(pid,corner)];
%      end     
%     end
%     Ncell = unique(Ncell);
%       
%     Number_cell(n) = length(Ncell);
    
%     % If out off boundary
%     if Number_cell(n) == 0
%         for i = 1:Np
%          pid = npoints{n}(1,i);
%          corner = npoints{n}(2,i);
%          id_n = find(CONNECT_MLS{pid,corner}==n);
% 
%          % Compute shape function from MLS
%          N_MLS{pid,corner}(id_n) = 0;
%         end
%     
%         % If the node in the object boundary
%         elseif Number_cell(n)<4 && Number_cell(n) > 0
%         color_node(n) = 2;
%         node_boundary = [node_boundary n];
%         
%         % If the node in the object boundary
%         N_total = 0;
%         for i = 1:Np
%          pid = npoints{n}(1,i);
%          corner = npoints{n}(2,i);
%          id_n = find(CONNECT_MLS{pid,corner}==n);
%          
% %          N_MLS{pid,corner}(id_n) = linearshape(x_data{pid,corner},LOC(n,:),le(1),le(2));
%          N_MLS{pid,corner}(id_n) = Cubic_Bspline(x_data{pid,corner},LOC(n,:),le(1),le(2));
%          
%          N_total = N_total + N_MLS{pid,corner}(id_n);
%         end
%         
%          for i = 1:Np
%          pid = npoints{n}(1,i);
%          corner = npoints{n}(2,i);
%          id_n = find(CONNECT_MLS{pid,corner}==n);
%              N_MLS{pid,corner}(id_n) = N_MLS{pid,corner}(id_n)/N_total;
%          end  
         
%         % If the node in the object body
%          elseif Number_cell(n) == 4
        color_node(n) = 1;
        
             cv = [];      
     won = ones(1,Np);
     wn = zeros(1,Np);
     
     for i = 1:Np
         pid = npoints{n}(1,i);
         corner = npoints{n}(2,i);
         
         cv = [cv x_data{pid,corner}];       
%          Nn_local = linearshape(x_data{pid,corner},LOC(n,:),le(1),le(2));
         Nn_local = Cubic_Bspline(x_data{pid,corner},LOC(n,:),le(1),le(2));
         
          wn(i) = Nn_local;
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


% Compute condition number
nb = 3;
j = 1;
while j<=nb
    if gama(j)==0
        gama(j)=[];
        nb = nb-1;
    end
    j = j+1;
end

k(n) = max(gama)/min(gama);
% if Np == 1
% %     phi(i)=wn(i);  
%     node_boundary = [node_boundary n];
%     
%     
%          for i = 1:Np
%          pid = npoints{n}(1,i);
%          corner = npoints{n}(2,i);
%          id_n = find(CONNECT_MLS{pid,corner}==n);
% 
%          % Compute shape function from MLS
%          N_MLS{pid,corner}(id_n) = 0;
%          end
%          
%              
% elseif k(n)>100
% %     for i=1:Np  
% %     phi(i)=wn(i)/gama(1);  
% %     end
%     node_boundary = [node_boundary n];
%         for i = 1:Np
%          pid = npoints{n}(1,i);
%          corner = npoints{n}(2,i);
%          id_n = find(CONNECT_MLS{pid,corner}==n);
% 
%          % Compute shape function from MLS
%          N_MLS{pid,corner}(id_n) = 0;
%          end
% else
    phi=zeros(1,Np);
    for i=1:Np
        cji=0;
        for j=1:3
            cji=cji+cc(j,i);
        end   
        phi(i)=wn(i)*cji;  
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for i = 1:Np
         pid = npoints{n}(1,i);
         corner = npoints{n}(2,i);
         id_n = find(CONNECT_MLS{pid,corner}==n);

         % Compute shape function from MLS
         N_MLS{pid,corner}(id_n) = phi(i);
         end
         
% end
 end