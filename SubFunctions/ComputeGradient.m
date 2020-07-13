function [DensityGradient,StressGradient,density_data,stress_data] = ComputeGradient...
    (spCount,nodeCount,cellCount,mspoints,CONNECT,N,dN,NODES,p_sp,V_sp,s_sp,x_sp,r1_sp,r2_sp)

% Compute the gradient term
ndensity            = zeros(nodeCount,1);
nvolume             = zeros(nodeCount,1);
nstress             = cell(nodeCount,1);

x_corner1 = x_sp - r1_sp - r2_sp;
x_corner2 = x_sp + r1_sp - r2_sp;
x_corner3 = x_sp + r1_sp + r2_sp;
x_corner4 = x_sp - r1_sp + r2_sp;


for n = 1:nodeCount
    nstress{n} = zeros(3,1);
end

for p=1:spCount
 for j=1:NODES(p)
     npid                           = CONNECT{p}(j);
     
          if N{p}(j)==1
         continue
          end
     
 % Density
 ndensity(npid)         = ndensity(npid) + p_sp(p)*V_sp(p)*N{p}(j);
 
 % Volume
 nvolume(npid)          = nvolume(npid,:) + V_sp(p)*N{p}(j);
 
 % Internal forces
 nstress{npid}          = nstress{npid} + (V_sp(p)*s_sp(p,:)*N{p}(:,j))';
 end 
 end

for n = 1:nodeCount
    if nvolume>0
        ndensity(n) = ndensity(n)/nvolume(n);
        nstress{n} = nstress{n}/nvolume(n);
    end
end

for c = 1:cellCount
    mpts = mspoints{c};
    
    for sp = 1:length(mpts)
        spid = mpts(sp);
        DensityGradient{spid} = zeros(1,2);
        StressGradient{spid} = zeros(3,2);
        
    for j=1:NODES(spid)
              if dN{spid}(j)==0
             continue
              end
                npid = CONNECT{spid}(j);
                DensityGradient{spid} = DensityGradient{spid} + (ndensity(npid)*dN{spid}(:,j)');
                StressGradient{spid} = StressGradient{spid} + (nstress{npid}*dN{spid}(:,j)');
    end 
    end      
end

density_data = cell(spCount,5);
stress_data = cell(spCount,5);

% Corner position and velocity
 for sp = 1:spCount         
    density_data{sp,1} = p_sp(sp);
    density_data{sp,2} = p_sp(sp) + DensityGradient{sp}*(x_corner1(sp,:)- x_sp(sp,:))';    
    density_data{sp,3} = p_sp(sp) + DensityGradient{sp}*(x_corner2(sp,:)- x_sp(sp,:))'; 
    density_data{sp,4} = p_sp(sp) + DensityGradient{sp}*(x_corner3(sp,:)- x_sp(sp,:))'; 
    density_data{sp,5} = p_sp(sp) + DensityGradient{sp}*(x_corner4(sp,:)- x_sp(sp,:))'; 
    
    stress_data{sp,1} = s_sp(sp,:);
    stress_data{sp,2} = s_sp(sp,:) + (StressGradient{sp}*(x_corner1(sp,:)- x_sp(sp,:))')';    
    stress_data{sp,3} = s_sp(sp,:) + (StressGradient{sp}*(x_corner2(sp,:)- x_sp(sp,:))')'; 
    stress_data{sp,4} = s_sp(sp,:) + (StressGradient{sp}*(x_corner3(sp,:)- x_sp(sp,:))')'; 
    stress_data{sp,5} = s_sp(sp,:) + (StressGradient{sp}*(x_corner4(sp,:)- x_sp(sp,:))')'; 
 end