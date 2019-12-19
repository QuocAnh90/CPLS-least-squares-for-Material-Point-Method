
% RMS Errors Vector 
MPM_displacement=zeros(6,1);MPM_velocity=zeros(6,1);MPM_stress=zeros(6,1);

BSMPM_displacement=zeros(6,1);BSMPM_velocity=zeros(6,1);BSMPM_stress=zeros(6,1);

MLS_BSMPM_displacement=zeros(6,1);MLS_BSMPM_velocity=zeros(6,1);MLS_BSMPM_stress=zeros(6,1);

% Run simulations
% RMS results vector store errors
% Vibration('interpolator',resolution)

% %% MPM interpolator
% results = Vibration('MPM',8);
% MPM_displacement(1) = results(1);MPM_velocity(1) = results(2);MPM_stress(1) = results(3);
% clear results
% results = Vibration('MPM',16);
% MPM_displacement(2) = results(1);MPM_velocity(2) = results(2);MPM_stress(2) = results(3);
% clear results
% results = Vibration('MPM',32);
% MPM_displacement(3) = results(1);MPM_velocity(3) = results(2);MPM_stress(3) = results(3);
% clear results
% results = Vibration('MPM',64);
% MPM_displacement(4) = results(1);MPM_velocity(4) = results(2);MPM_stress(4) = results(3);
% clear results
% 
% % BSMPM interpolator
% results = Vibration('BSMPM',8);
% BSMPM_displacement(1) = results(1);BSMPM_velocity(1) = results(2);BSMPM_stress(1) = results(3);
% clear results
% results = Vibration('BSMPM',16);
% BSMPM_displacement(2) = results(1);BSMPM_velocity(2) = results(2);BSMPM_stress(2) = results(3);
% clear results
% results = Vibration('BSMPM',32);
% BSMPM_displacement(3) = results(1);BSMPM_velocity(3) = results(2);BSMPM_stress(3) = results(3);
% clear results
% results = Vibration('BSMPM',64);
% BSMPM_displacement(4) = results(1);BSMPM_velocity(4) = results(2);BSMPM_stress(4) = results(3);
% clear results
% results = Vibration('BSMPM',128);
% BSMPM_displacement(5) = results(1);BSMPM_velocity(5) = results(2);BSMPM_stress(5) = results(3);
% clear results
% results = Vibration('BSMPM',256);
% BSMPM_displacement(6) = results(1);BSMPM_velocity(6) = results(2);BSMPM_stress(6) = results(3);
% clear results

%% MLS_BSMPM interpolator
results = Vibration('MLS_BSMPM',8);
MLS_BSMPM_displacement(1) = results(1);MLS_BSMPM_velocity(1) = results(2);MLS_BSMPM_stress(1) = results(3);
clear results
results = Vibration('MLS_BSMPM',16);
MLS_BSMPM_displacement(2) = results(1);MLS_BSMPM_velocity(2) = results(2);MLS_BSMPM_stress(2) = results(3);
clear results
results = Vibration('MLS_BSMPM',32);
MLS_BSMPM_displacement(3) = results(1);MLS_BSMPM_velocity(3) = results(2);MLS_BSMPM_stress(3) = results(3);
clear results
results = Vibration('MLS_BSMPM',64);
MLS_BSMPM_displacement(4) = results(1);MLS_BSMPM_velocity(4) = results(2);MLS_BSMPM_stress(4) = results(3);
clear results
% results = Vibration('MLS_BSMPM',128);
% MLS_BSMPM_displacement(5) = results(1);MLS_BSMPM_velocity(5) = results(2);MLS_BSMPM_stress(5) = results(3);
% clear results
% results = Vibration('MLS_BSMPM',256);
% MLS_BSMPM_displacement(6) = results(1);MLS_BSMPM_velocity(6) = results(2);MLS_BSMPM_stress(6) = results(3);
% clear results

%% Plot results
xconverge1 = [0.0316 0.1 0.1 0.0316];
yconverge1 = [0.0000002 0.0000002 0.000002 0.0000002];
% yconverge1 = [0.0002 0.0002 0.002 0.0002];
text2 ='2';
text1 = '1';
sz = 100;

dcell = [0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625];
figure 
loglog(dcell,BSMPM_displacement,'g',dcell,MLS_BSMPM_displacement,'--b','linewidth',3)
hold on
plot(xconverge1,yconverge1,'k','linewidth',2)
hold on
scatter(dcell,MPM_displacement,sz,'r','filled')
hold on
scatter(dcell,BSMPM_displacement,sz,'g','filled')
hold on
scatter(dcell,MLS_BSMPM_displacement,sz,'b','filled')
text(0.11,0.0000007,text2,'FontSize',18)
text(0.055,0.00000015,text1,'FontSize',18)
% text(0.11,0.0007,text2,'FontSize',18)
% text(0.055,0.00015,text1,'FontSize',18)
title('displacement converegene')
legend('BSMPM','BSMPM with MLS')
xlabel('cell size (m)')
ylabel('RMS')
% axis([0.01 0.2 0.0001 0.1])
set(gca,'fontsize', 18)

xconverge2 = [0.0316 0.1 0.1 0.0316];
yconverge2 = [0.2 0.2 2 0.2];
figure 
loglog(dcell,MPM_velocity,'r',dcell,BSMPM_velocity,'g',dcell,MLS_BSMPM_velocity,'--b','linewidth',3)
hold on
% plot(xconverge2,yconverge2,'k','linewidth',2)
% hold on
scatter(dcell,MPM_velocity,sz,'r','filled')
hold on
scatter(dcell,BSMPM_velocity,sz,'g','filled')
hold on
scatter(dcell,MLS_BSMPM_velocity,sz,'b','filled')
% text(0.11,0.7,text2,'FontSize',18)
% text(0.055,0.15,text1,'FontSize',18)
% title('velocity converegene')
legend('MPM','BSMPM','BSMPM with MLS')
xlabel('cell size (m)')
ylabel('RMS')
% axis([0.01 0.2 0.1 100])

xconverge3 = [0.0316 0.1 0.1 0.0316];
yconverge3 = [20000 20000 200000 20000];
figure 
loglog(dcell,MPM_stress,'r',dcell,BSMPM_stress,'g',dcell,MLS_BSMPM_stress,'--b','linewidth',3)
hold on
% plot(xconverge3,yconverge3,'k','linewidth',2)
% hold on
scatter(dcell,MPM_stress,sz,'r','filled')
hold on
scatter(dcell,BSMPM_stress,sz,'g','filled')
hold on
scatter(dcell,MLS_BSMPM_stress,sz,'b','filled')
% text(0.11,70000,text2,'FontSize',18)
% text(0.055,15000,text1,'FontSize',18)
% title('stress converegene')
legend('MPM','BSMPM','BSMPM with MLS')
xlabel('cell size (m)')
ylabel('RMS')
% axis([0.01 0.2 10000 10000000])