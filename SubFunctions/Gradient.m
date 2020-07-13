function [velocity_data,x_data] = Gradient(spCount,x_sp,r1_sp,r2_sp,v_ssp,L_sp)


%% Update velocity field at the corner
velocity_data = cell(spCount,5);
x_data = cell(spCount,5);

x_corner1 = x_sp - r1_sp - r2_sp;
x_corner2 = x_sp + r1_sp - r2_sp;
x_corner3 = x_sp + r1_sp + r2_sp;
x_corner4 = x_sp - r1_sp + r2_sp;

% Corner position and velocity
 for sp = 1:spCount
    x_data{sp,1} = x_sp(sp,:)';             
    velocity_data{sp,1} = v_ssp(sp,:)';
    x_data{sp,2} = x_corner1(sp,:)';        
    velocity_data{sp,2} = v_ssp(sp,:)' + L_sp{sp}*(x_corner1(sp,:)- x_sp(sp,:))';    
    x_data{sp,3} = x_corner2(sp,:)';        
    velocity_data{sp,3} = v_ssp(sp,:)' + L_sp{sp}*(x_corner2(sp,:)- x_sp(sp,:))'; 
    x_data{sp,4} = x_corner3(sp,:)';        
    velocity_data{sp,4} = v_ssp(sp,:)' + L_sp{sp}*(x_corner3(sp,:)- x_sp(sp,:))'; 
    x_data{sp,5} = x_corner4(sp,:)';        
    velocity_data{sp,5} = v_ssp(sp,:)' + L_sp{sp}*(x_corner4(sp,:)- x_sp(sp,:))'; 
 end