function [N_g] = CellMLS...
    (spElems_boundary,spCount,cellCount,cpoints,x_data,LOCC,le,Elems_CPDI,CONNECT_MLS_G)

 %% Quadratic-Bspline for mapping from particle to centroids using IMLS
  N_g         = cell(spCount,5);

 for c = 1:cellCount
     if isempty(cpoints{c})==1
         continue
     end
      Np = length(cpoints{c}(1,:));
      
     if ismember(c,spElems_boundary)==1
         
     for i = 1:Np
         pid = cpoints{c}(1,i);
         corner = cpoints{c}(2,i);
         spElems_test = Elems_CPDI(pid,corner);
         id_g = find(CONNECT_MLS_G{pid,corner}==c);        
             N_g{pid,corner}(id_g) = 0;
     end
     
     elseif ismember(c,spElems_boundary)==0
     cv = [];
     won = ones(1,Np);
     w = zeros(1,Np);
     
     for i = 1:Np
         pid = cpoints{c}(1,i);
         corner = cpoints{c}(2,i);
         cv = [cv x_data{pid,corner}];
         
         N_local = Quadratic_Bspline(x_data{pid,corner},LOCC(c,:),le(1),le(2));   
         w(i) = N_local;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Compute MLS shape function using cell boundaries     
     for i = 1:Np
         pid = cpoints{c}(1,i);
         corner = cpoints{c}(2,i);
         spElems_test = Elems_CPDI(pid,corner);
         id_g = find(CONNECT_MLS_G{pid,corner}==c);        
             N_g{pid,corner}(id_g) = phi(i);
     end
     
     end
  end