clear all
close all

% RMS Errors Vector 
CPDI_displacement=zeros(7,1);CPDI_velocity=zeros(7,1); CPDI_stress=zeros(7,1);

IMLS_CPDI_displacement=zeros(7,1);IMLS_CPDI_velocity=zeros(7,1); IMLS_CPDI_stress=zeros(7,1);

% Run simulations
% RMS results vector store errors
%% CPDI interpolator
% % results = Vortex('CPDI',1); 
% % CPDI_displacement(1) = results(1);CPDI_velocity(1) = results(2); CPDI_stress(1) = results(3);
% % clear results
results = Vortex('CPDI',2);
CPDI_displacement(2) = results(1);CPDI_velocity(2) = results(2); CPDI_stress(2) = results(3);
clear results
results = Vortex('CPDI',4);
CPDI_displacement(3) = results(1);CPDI_velocity(3) = results(2); CPDI_stress(3) = results(3);
clear results
results = Vortex('CPDI',8);
CPDI_displacement(4) = results(1);CPDI_velocity(4) = results(2); CPDI_stress(4) = results(3);
clear results
results = Vortex('CPDI',16);
CPDI_displacement(5) = results(1);CPDI_velocity(5) = results(2);CPDI_stress(5) = results(3);
clear results
% results = Vortex('CPDI',32);
% CPDI_displacement(6) = results(1);CPDI_velocity(6) = results(2);CPDI_stress(6) = results(3);
% clear results
% results = Vortex('CPDI',64);
% CPDI_displacement(7) = results(1);CPDI_velocity(7) = results(2);CPDI_stress(5) = results(3);
% clear results

%% IMLS_CPDI interpolator
% results = Vortex('IMLS_CPDI',1);
% IMLS_CPDI_displacement(1) = results(1);IMLS_CPDI_velocity(1) = results(2); IMLS_CPDI_stress(1) = results(3);
% % clear results
results = Vortex('IMLS_CPDI',2);
IMLS_CPDI_displacement(2) = results(1);IMLS_CPDI_velocity(2) = results(2); IMLS_CPDI_stress(2) = results(3);
clear results
results = Vortex('IMLS_CPDI',4);
IMLS_CPDI_displacement(3) = results(1);IMLS_CPDI_velocity(3) = results(2); IMLS_CPDI_stress(3) = results(3);
clear results
results = Vortex('IMLS_CPDI',8);
IMLS_CPDI_displacement(4) = results(1);IMLS_CPDI_velocity(4) = results(2); IMLS_CPDI_stress(4) = results(3);
clear results
results = Vortex('IMLS_CPDI',16);
IMLS_CPDI_displacement(5) = results(1);IMLS_CPDI_velocity(5) = results(2); IMLS_CPDI_stress(5) = results(3);
clear results
% results = Vortex('IMLS_CPDI',32);
% IMLS_CPDI_displacement(6) = results(1);IMLS_CPDI_velocity(6) = results(2); IMLS_CPDI_stress(6) = results(3);
% clear results
% % results = Vortex('IMLS_CPDI',64);
% % IMLS_CPDI_displacement(7) = results(1);IMLS_CPDI_velocity(7) = results(2);% IMLS_CPDI_stress(4) = results(3);
% % clear results

%% Plot results
dcell = [0.25 0.125 0.0625 0.03125 0.01563 0.01563/2 0.01563/4];
figure 
loglog(dcell,CPDI_displacement,'--dg',dcell,IMLS_CPDI_displacement,'-sb','linewidth',3)
hold on
x1 = [0.02 0.2]; y1 = [0.002 0.2];
plot(x1,y1,'--k','linewidth',2);
hold on
% scatter(dcell,CPDI_displacement,100,'g','filled')
% hold on
% scatter(dcell,IMLS_CPDI_displacement,100,'b','filled')

% title('displacement converegene')
legend('CPDI','CPLS','2nd order','Location','southeast')
xlabel('cell size (m)')
ylabel('RMS')
set(gca,'fontsize', 15)

figure 
loglog(dcell,CPDI_velocity,'--dg',dcell,IMLS_CPDI_velocity,'-sb','linewidth',3)
hold on
x2 = [0.03 0.3]; y2 = [0.02 2];
plot(x2,y2,'--k','linewidth',2);
hold on
% scatter(dcell,CPDI_velocity,100,'g','filled')
% hold on
% scatter(dcell,IMLS_CPDI_velocity,100,'b','filled')

% title('velocity converegene')
legend('CPDI','CPLS','2nd order','Location','southeast')
xlabel('cell size (m)')
ylabel('RMS')
set(gca,'fontsize', 15)

figure 
loglog(dcell,CPDI_stress,'--dg',dcell,IMLS_CPDI_stress,'-sb','linewidth',3)
hold on
x3 = [0.03 0.3]; y3 = [20 2000];
plot(x3,y3,'--k','linewidth',2);
hold on
% scatter(dcell,CPDI_stress,100,'g','filled')
% hold on
% scatter(dcell,IMLS_CPDI_stress,100,'b','filled')

% title('stress converegene')
legend('CPDI','CPLS','2nd order','Location','southeast')
xlabel('cell size (m)')
ylabel('RMS')
set(gca,'fontsize', 15)

