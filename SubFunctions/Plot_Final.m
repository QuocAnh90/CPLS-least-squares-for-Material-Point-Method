function StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,color_node,v_ssp,spCount,spElems_boundary,cellCount,nodeCount,r1_sp,r2_sp,p_sc,LOCC)

 x_corner1              = zeros(spCount,2);
 x_corner2              = zeros(spCount,2);
 x_corner3              = zeros(spCount,2);
 x_corner4              = zeros(spCount,2);
 color_cell             = zeros(cellCount,1);
 
displacement = zeros(spCount,1);
velocity = zeros(spCount,1);
for sp=1:spCount
    displacement(sp) = sqrt(d_sp(sp,1)^2+d_sp(sp,2)^2);
    velocity(sp) = sqrt(v_ssp(sp,1)^2+v_ssp(sp,2)^2);
end

 for sp=1:spCount
 x_corner1(sp,:) = x_sp(sp,:) - r1_sp(sp,:) - r2_sp(sp,:);
 x_corner2(sp,:) = x_sp(sp,:) + r1_sp(sp,:) - r2_sp(sp,:);
 x_corner3(sp,:) = x_sp(sp,:) + r1_sp(sp,:) + r2_sp(sp,:);
 x_corner4(sp,:) = x_sp(sp,:) - r1_sp(sp,:) + r2_sp(sp,:);
 end

StressProfile1=figure;
% figure
% set(StressProfile1, 'visible','off');
sz = 200*le(1);
color = velocity;

%% Cell color
 for c = 1:cellCount
    if p_sc(c)>0
        color_cell(c) = max(color);
    end
    
    if ismember(c,spElems_boundary)==1
        color_cell(c) = max(color)/2;
    end
end

% Remove zero color node
x_cell = LOCC;
cellCount1 = cellCount;

j=1;
while j < cellCount1
    if color_cell(j) == 0
        color_cell(j) = [];
        x_cell(j,:) = [];
        cellCount1 = cellCount1-1;
        continue
    end
    j = j+1;
end
 
% Node color
for n = 1:length(color_node)
    if color_node(n) ==2
        color_node(n) = 1.2;
    elseif color_node(n) ==1
        color_node(n) = 3;
    else
        color_node(n) = 0;
    end
end


% Remove zero color node
x_node = LOC;
nodeCount1 = nodeCount;
n=1;
while n < nodeCount1
    if color_node(n) == 0
        color_node(n) = [];
        x_node(n,:) = [];
        nodeCount1 = nodeCount1-1;
        continue
    end
    n = n+1;
end

%% Start ploting

%% Plot particle
scatter(x_sp(:,1),x_sp(:,2),sz,color,'filled');
hold on

%% Plot corner point
% scatter(x_corner1(:,1),x_corner1(:,2));
% hold on
% scatter(x_corner2(:,1),x_corner2(:,2));
% hold on
% scatter(x_corner3(:,1),x_corner3(:,2));
% hold on
% scatter(x_corner4(:,1),x_corner4(:,2));
% hold on

%% Plot particle domain
for sp=1:spCount
% %     Qxx{sp} = [x_sp(sp,1) x_sp(sp,1)+r1_sp(sp,1)];
% %     Qxy{sp} = [x_sp(sp,2) x_sp(sp,2)+r1_sp(sp,2)];
% %     Qyx{sp} = [x_sp(sp,1) x_sp(sp,1)+r2_sp(sp,1)];
% %     Qyy{sp} = [x_sp(sp,2) x_sp(sp,2)+r2_sp(sp,2)];
% %     plot (Qxx{sp},Qxy{sp},'r')
% %     hold on 
% %     plot (Qyx{sp},Qyy{sp},'r')
% %     hold on
%     
%     x_cor = [x_corner1(sp,1) x_corner2(sp,1) x_corner3(sp,1) x_corner4(sp,1) x_corner1(sp,1)];
%     y_cor = [x_corner1(sp,2) x_corner2(sp,2) x_corner3(sp,2) x_corner4(sp,2) x_corner1(sp,2)];
%     plot(x_cor,y_cor,'b')
%     hold on
end

% %% Determine the cell boundary of object
%     Ncorner_boundary = length(boundary_corner(1,:));
% 
%  for i = 1:Ncorner_boundary
%      index_corner = boundary_corner(1:2,i);
%      particle_index = index_corner(1);
%      corner_index = index_corner(2);
%      
%      scatter(x_data{particle_index,corner_index}(1),x_data{particle_index,corner_index}(2),10,'r','filled');
%      hold on
%  end
 
% Plot node color
% scatter(x_node(:,1),x_node(:,2),sz,color_node,'filled');
% hold on

% Plot cell color
% scatter(x_cell(:,1),x_cell(:,2),sz,color_cell,'filled');
% hold on
grid on

%% Fornat
axis([0,max(LOC(:,1)),0,max(LOC(:,2))]);
% axis([1.25,3.75,1.25,3.75]);
% axis([0.9,1.5,1.5,3]);
set(gca,'xtick',[0:le(1):max(LOC(:,1))]);
set(gca,'ytick',[0:le(2):max(LOC(:,2))]);
h=colorbar;
colormap(jet(256))
% set(h, 'ytick', [0:0.2:1.2]);
    caxis([0 0.1]);