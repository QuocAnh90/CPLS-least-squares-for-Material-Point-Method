clear all
close all

% RMS Errors Vector 
MPM_displacement=zeros(6,1);MPM_velocity=zeros(6,1);MPM_stress=zeros(6,1);

CPDI_displacement=zeros(6,1);CPDI_velocity=zeros(6,1);CPDI_stress=zeros(6,1);
% 
IMLS_CPDI_displacement=zeros(6,1);IMLS_CPDI_velocity=zeros(6,1);IMLS_CPDI_stress=zeros(6,1);


%% Run simulations
% RMS results vector store errors
% Vibration('interpolator',resolution)

% MPM interpolator
results = Vibration2D('MPM',8);
MPM_displacement(1) = results(1);MPM_velocity(1) = results(2);MPM_stress(1) = results(3);
clear results
results = Vibration2D('MPM',16);
MPM_displacement(2) = results(1);MPM_velocity(2) = results(2);MPM_stress(2) = results(3);
% clear results
% results = Vibration2D('MPM',32);
% MPM_displacement(3) = results(1);MPM_velocity(3) = results(2);MPM_stress(3) = results(3);
% clear results
% results = Vibration2D('MPM',64);
% MPM_displacement(4) = results(1);MPM_velocity(4) = results(2);MPM_stress(4) = results(3);
% clear results
% results = Vibration2D('MPM',128);
% MPM_displacement(5) = results(1);MPM_velocity(5) = results(2);MPM_stress(5) = results(3);
% clear results

%% CPDI interpolator
results = Vibration2D('CPDI',8);
CPDI_displacement(1) = results(1);CPDI_velocity(1) = results(2);CPDI_stress(1) = results(3);
clear results
results = Vibration2D('CPDI',16);
CPDI_displacement(2) = results(1);CPDI_velocity(2) = results(2);CPDI_stress(2) = results(3);
clear results
results = Vibration2D('CPDI',32);
CPDI_displacement(3) = results(1);CPDI_velocity(3) = results(2);CPDI_stress(3) = results(3);
clear results
results = Vibration2D('CPDI',64);
CPDI_displacement(4) = results(1);CPDI_velocity(4) = results(2);CPDI_stress(4) = results(3);
clear results
results = Vibration2D('CPDI',128);
CPDI_displacement(5) = results(1);CPDI_velocity(5) = results(2);CPDI_stress(5) = results(3);
clear results
% results = Vibration2D('CPDI',256);
% CPDI_displacement(6) = results(1);CPDI_velocity(6) = results(2);CPDI_stress(6) = results(3);
% clear results

%% IMLS_CPDI interpolator
results = Vibration2D('IMLS_CPDI',8);
IMLS_CPDI_displacement(1) = results(1);IMLS_CPDI_velocity(1) = results(2);IMLS_CPDI_stress(1) = results(3);
clear results
results = Vibration2D('IMLS_CPDI',16);
IMLS_CPDI_displacement(2) = results(1);IMLS_CPDI_velocity(2) = results(2);IMLS_CPDI_stress(2) = results(3);
clear results
results = Vibration2D('IMLS_CPDI',32);
IMLS_CPDI_displacement(3) = results(1);IMLS_CPDI_velocity(3) = results(2);IMLS_CPDI_stress(3) = results(3);
clear results
results = Vibration2D('IMLS_CPDI',64);
IMLS_CPDI_displacement(4) = results(1);IMLS_CPDI_velocity(4) = results(2);IMLS_CPDI_stress(4) = results(3);
clear results
results = Vibration2D('IMLS_CPDI',128);
IMLS_CPDI_displacement(5) = results(1);IMLS_CPDI_velocity(5) = results(2);IMLS_CPDI_stress(5) = results(3);
clear results
% results = Vibration2D('IMLS_CPDI',256);
% IMLS_CPDI_displacement(6) = results(1);IMLS_CPDI_velocity(6) = results(2);IMLS_CPDI_stress(6) = results(3);
% clear results

%% Plot results
dcell = [0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625];
figure 
loglog(dcell(1:4),MPM_displacement(1:4),'r',dcell,CPDI_displacement,'--g',dcell,IMLS_CPDI_displacement,'--b','linewidth',3)
hold on
scatter(dcell(1:4),MPM_displacement(1:4),100,'r','filled')
hold on
scatter(dcell,CPDI_displacement,100,'d','g','filled')
hold on
scatter(dcell,IMLS_CPDI_displacement,100,'s','b','filled')
hold on
title('displacement converegene')
legend('MPM','CPDI','CPDI with IMLS','Location','southeast')
xlabel('log(le)')
ylabel('RMS')

figure 
loglog(dcell(1:4),MPM_velocity(1:4),'r',dcell,CPDI_velocity,'--g',dcell,IMLS_CPDI_velocity,'--b','linewidth',3)
hold on
scatter(dcell(1:4),MPM_velocity(1:4),100,'r','filled')
hold on
scatter(dcell,CPDI_velocity,100,'d','g','filled')
hold on
scatter(dcell,IMLS_CPDI_velocity,100,'s','b','filled')
hold on
title('velocity converegene')
legend('MPM','CPDI','CPDI with IMLS','Location','southeast')
xlabel('log(le)')
ylabel('RMS')

figure 
loglog(dcell(1:4),MPM_stress(1:4),'r',dcell,CPDI_stress,'--g',dcell,IMLS_CPDI_stress,'--b','linewidth',3)
hold on
scatter(dcell(1:4),MPM_stress(1:4),100,'r','filled')
hold on
scatter(dcell,CPDI_stress,100,'d','g','filled')
hold on
scatter(dcell,IMLS_CPDI_stress,100,'s','b','filled')
hold on
title('stress converegene')
legend('MPM','CPDI','CPDI with IMLS','Location','southeast')
xlabel('log(le)')
ylabel('RMS')