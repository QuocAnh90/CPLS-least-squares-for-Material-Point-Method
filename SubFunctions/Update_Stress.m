function [F_sp,V_sp,s_sp,p_sp,L_sp] = Update_Stress(CModel,CModel_parameter,...
    NODES,dt,cellCount,mspoints,CONNECT,nvelo_si,dN,...
    F_sp,V_spo,m_sp,s_sp,p_sp,V_sp)

% Calculate stress for solid phase
for c = 1:cellCount
    mpts = mspoints{c};
    
    for sp = 1:length(mpts)
        spid = mpts(sp);
        L_sp{spid} = zeros(2,2);
    
for j=1:NODES(spid)
%           if dN{spid}(j)==0
%          continue
%           end
            npid = CONNECT{spid}(j);
            L_sp{spid} = L_sp{spid} + (nvelo_si(npid,:)'*dN{spid}(:,j)');
end          
        dESP = (L_sp{spid} + L_sp{spid}')/2*dt;
        
%         % Cycle for deformation gradient approximation
%         F = F_sp{spid};
%         cycle = 1000;
%         dt_cycle = dt/cycle;
%         
%         dF = (eye(2,2)+L_sp{spid}*dt_cycle);
%         
%         for i = 1:cycle
%             F = dF * F; 
%         end
%         F_sp{spid} = F;
        
        F_sp{spid} = (eye(2,2)+L_sp{spid}*dt)*F_sp{spid};    
        
        J = det(F_sp{spid});
        V_sp(spid)=V_spo(spid)*J;   

        switch CModel
            case 'Neo_Hookean_Elastic'
                [s_sp(spid,:)]=Neo_Hookean_elastic(CModel_parameter,F_sp{spid},J);
            case 'Linear_Elastic'
                [s_sp(spid,:)]=Linear_elastic(CModel_parameter,dESP,s_sp(spid,:));             
            case 'Water'
                [s_sp(spid,:)]=Water(CModel_parameter,(L_sp{spid} + L_sp{spid}')/2,J);
        end
        p_sp(spid) = m_sp(spid)/V_sp(spid);
    end      
end