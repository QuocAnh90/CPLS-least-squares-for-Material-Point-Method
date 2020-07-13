% clear all
% close all

% RMS Errors Vector 
% MPM_displacement=zeros(7,1);MPM_velocity=zeros(7,1);MPM_stress=zeros(7,1);
% MPM_displacement4=zeros(7,1);MPM_velocity4=zeros(7,1);MPM_stress4=zeros(7,1);
% MPM_displacement9=zeros(7,1);MPM_velocity9=zeros(7,1);MPM_stress9=zeros(7,1);


CPDI_displacement=zeros(7,1);CPDI_velocity=zeros(7,1);CPDI_stress=zeros(7,1);
CPDI_displacement4=zeros(7,1);CPDI_velocity4=zeros(7,1);CPDI_stress4=zeros(7,1);
CPDI_displacement9=zeros(7,1);CPDI_velocity9=zeros(7,1);CPDI_stress9=zeros(7,1);

% IMLS_MPM_displacement=zeros(7,1);IMLS_MPM_velocity=zeros(7,1);IMLS_MPM_stress=zeros(7,1);


% gradientCPDI_displacement=zeros(7,1);gradientCPDI_velocity=zeros(7,1);gradientCPDI_stress=zeros(7,1);

% IMLS_CPDI_displacement=zeros(7,1);IMLS_CPDI_velocity=zeros(7,1);IMLS_CPDI_stress=zeros(7,1);
% IMLS_CPDI_displacement4=zeros(7,1);IMLS_CPDI_velocity4=zeros(7,1);IMLS_CPDI_stress4=zeros(7,1);
% IMLS_CPDI_displacement9=zeros(7,1);IMLS_CPDI_velocity9=zeros(7,1);IMLS_CPDI_stress9=zeros(7,1);



%% Run simulations
% RMS results vector store errors
% Vibration('interpolator',resolution)

% % MPM interpolator
% results = Vibration('MPM',8,1);
% MPM_displacement(1) = results(1);MPM_velocity(1) = results(2);MPM_stress(1) = results(3);
% clear results
% results = Vibration('MPM',16,1);
% MPM_displacement(2) = results(1);MPM_velocity(2) = results(2);MPM_stress(2) = results(3);
% clear results
% results = Vibration('MPM',32,1);
% MPM_displacement(3) = results(1);MPM_velocity(3) = results(2);MPM_stress(3) = results(3);
% clear results
% results = Vibration('MPM',64,1);
% MPM_displacement(4) = results(1);MPM_velocity(4) = results(2);MPM_stress(4) = results(3);
% clear results
% results = Vibration('MPM',128,1);
% MPM_displacement(5) = results(1);MPM_velocity(5) = results(2);MPM_stress(5) = results(3);
% clear results
% results = Vibration('MPM',256,1);
% MPM_displacement(6) = results(1);MPM_velocity(6) = results(2);MPM_stress(6) = results(3);
% clear results
% results = Vibration('MPM',512,1);
% MPM_displacement(7) = results(1);MPM_velocity(7) = results(2);MPM_stress(7) = results(3);
% clear results

% results = Vibration('MPM',8,4);
% MPM_displacement4(1) = results(1);MPM_velocity4(1) = results(2);MPM_stress4(1) = results(3);
% clear results
% results = Vibration('MPM',16,4);
% MPM_displacement4(2) = results(1);MPM_velocity4(2) = results(2);MPM_stress4(2) = results(3);
% clear results
% results = Vibration('MPM',32,4);
% MPM_displacement4(3) = results(1);MPM_velocity4(3) = results(2);MPM_stress4(3) = results(3);
% clear results
% results = Vibration('MPM',64,4);
% MPM_displacement4(4) = results(1);MPM_velocity4(4) = results(2);MPM_stress4(4) = results(3);
% clear results
% results = Vibration('MPM',128,4);
% MPM_displacement4(5) = results(1);MPM_velocity4(5) = results(2);MPM_stress4(5) = results(3);
% clear results
% results = Vibration('MPM',256,4);
% MPM_displacement4(6) = results(1);MPM_velocity4(6) = results(2);MPM_stress4(6) = results(3);
% clear results
% results = Vibration('MPM',512,4);
% MPM_displacement4(7) = results(1);MPM_velocity4(7) = results(2);MPM_stress4(7) = results(3);
% clear results

% results = Vibration('MPM',8,9);
% MPM_displacement9(1) = results(1);MPM_velocity9(1) = results(2);MPM_stress9(1) = results(3);
% clear results
% results = Vibration('MPM',16,9);
% MPM_displacement9(2) = results(1);MPM_velocity9(2) = results(2);MPM_stress9(2) = results(3);
% clear results
% results = Vibration('MPM',32,9);
% MPM_displacement9(3) = results(1);MPM_velocity9(3) = results(2);MPM_stress9(3) = results(3);
% clear results
% results = Vibration('MPM',64,9);
% MPM_displacement9(4) = results(1);MPM_velocity9(4) = results(2);MPM_stress9(4) = results(3);
% clear results
% results = Vibration('MPM',128,9);
% MPM_displacement9(5) = results(1);MPM_velocity9(5) = results(2);MPM_stress9(5) = results(3);
% clear results
% results = Vibration('MPM',256,9);
% MPM_displacement9(6) = results(1);MPM_velocity9(6) = results(2);MPM_stress9(6) = results(3);
% clear results
% results = Vibration('MPM',512,9);
% MPM_displacement9(7) = results(1);MPM_velocity9(7) = results(2);MPM_stress9(7) = results(3);
% clear results

% % % IMLS_MPM interpolator
% % results = Vibration('IMLS_MPM',8);
% % IMLS_MPM_displacement(1) = results(1);IMLS_MPM_velocity(1) = results(2);IMLS_MPM_stress(1) = results(3);
% % clear results
% % results = Vibration('IMLS_MPM',16);
% % IMLS_MPM_displacement(2) = results(1);IMLS_MPM_velocity(2) = results(2);IMLS_MPM_stress(2) = results(3);
% % clear results
% % results = Vibration('IMLS_MPM',32);
% % IMLS_MPM_displacement(3) = results(1);IMLS_MPM_velocity(3) = results(2);IMLS_MPM_stress(3) = results(3);
% % clear results
% % results = Vibration('IMLS_MPM',64);
% % IMLS_MPM_displacement(4) = results(1);IMLS_MPM_velocity(4) = results(2);IMLS_MPM_stress(4) = results(3);
% % clear results
% % results = Vibration('IMLS_MPM',128);
% % IMLS_MPM_displacement(5) = results(1);IMLS_MPM_velocity(5) = results(2);IMLS_MPM_stress(5) = results(3);
% % clear results
% % results = Vibration('IMLS_MPM',256);
% % IMLS_MPM_displacement(6) = results(1);IMLS_MPM_velocity(6) = results(2);IMLS_MPM_stress(6) = results(3);
% % clear results
% % results = Vibration('IMLS_MPM',512);
% % IMLS_MPM_displacement(7) = results(1);IMLS_MPM_velocity(7) = results(2);IMLS_MPM_stress(7) = results(3);
% % clear results
% 
% % % gradient CPDI interpolator
% % results = Vibration('CPDI',8);
% % gradientCPDI_displacement(1) = results(1);gradientCPDI_velocity(1) = results(2);gradientCPDI_stress(1) = results(3);
% % clear results
% % results = Vibration('CPDI',16);
% % gradientCPDI_displacement(2) = results(1);gradientCPDI_velocity(2) = results(2);gradientCPDI_stress(2) = results(3);
% % clear results
% % results = Vibration('CPDI',32);
% % gradientCPDI_displacement(3) = results(1);gradientCPDI_velocity(3) = results(2);gradientCPDI_stress(3) = results(3);
% % clear results
% % results = Vibration('CPDI',64);
% % gradientCPDI_displacement(4) = results(1);gradientCPDI_velocity(4) = results(2);gradientCPDI_stress(4) = results(3);
% % clear results
% % results = Vibration('CPDI',128);
% % gradientCPDI_displacement(5) = results(1);gradientCPDI_velocity(5) = results(2);gradientCPDI_stress(5) = results(3);
% % clear results
% % results = Vibration('CPDI',256);
% % gradientCPDI_displacement(6) = results(1);gradientCPDI_velocity(6) = results(2);gradientCPDI_stress(6) = results(3);
% % clear results
% % results = Vibration('CPDI',512);
% % gradientCPDI_displacement(7) = results(1);gradientCPDI_velocity(7) = results(2);gradientCPDI_stress(7) = results(3);
% % clear results
% 
% 
% CPDI interpolator
results = Vibration('CPDI',8,1);
CPDI_displacement(1) = results(1);CPDI_velocity(1) = results(2);CPDI_stress(1) = results(3);
clear results
results = Vibration('CPDI',16,1);
CPDI_displacement(2) = results(1);CPDI_velocity(2) = results(2);CPDI_stress(2) = results(3);
clear results
results = Vibration('CPDI',32,1);
CPDI_displacement(3) = results(1);CPDI_velocity(3) = results(2);CPDI_stress(3) = results(3);
clear results
results = Vibration('CPDI',64,1);
CPDI_displacement(4) = results(1);CPDI_velocity(4) = results(2);CPDI_stress(4) = results(3);
clear results
results = Vibration('CPDI',128,1);
CPDI_displacement(5) = results(1);CPDI_velocity(5) = results(2);CPDI_stress(5) = results(3);
clear results
results = Vibration('CPDI',256,1);
CPDI_displacement(6) = results(1);CPDI_velocity(6) = results(2);CPDI_stress(6) = results(3);
clear results
results = Vibration('CPDI',512,1);
CPDI_displacement(7) = results(1);CPDI_velocity(7) = results(2);CPDI_stress(7) = results(3);
clear results

results = Vibration('CPDI',8,4);
CPDI_displacement4(1) = results(1);CPDI_velocity4(1) = results(2);CPDI_stress4(1) = results(3);
clear results
results = Vibration('CPDI',16,4);
CPDI_displacement4(2) = results(1);CPDI_velocity4(2) = results(2);CPDI_stress4(2) = results(3);
clear results
results = Vibration('CPDI',32,4);
CPDI_displacement4(3) = results(1);CPDI_velocity4(3) = results(2);CPDI_stress4(3) = results(3);
clear results
results = Vibration('CPDI',64,4);
CPDI_displacement4(4) = results(1);CPDI_velocity4(4) = results(2);CPDI_stress4(4) = results(3);
clear results
results = Vibration('CPDI',128,4);
CPDI_displacement4(5) = results(1);CPDI_velocity4(5) = results(2);CPDI_stress4(5) = results(3);
clear results
results = Vibration('CPDI',256,4);
CPDI_displacement4(6) = results(1);CPDI_velocity4(6) = results(2);CPDI_stress4(6) = results(3);
clear results
results = Vibration('CPDI',512,4);
CPDI_displacement4(7) = results(1);CPDI_velocity4(7) = results(2);CPDI_stress4(7) = results(3);
clear results

results = Vibration('CPDI',8,9);
CPDI_displacement9(1) = results(1);CPDI_velocity9(1) = results(2);CPDI_stress9(1) = results(3);
clear results
results = Vibration('CPDI',16,9);
CPDI_displacement9(2) = results(1);CPDI_velocity9(2) = results(2);CPDI_stress9(2) = results(3);
clear results
results = Vibration('CPDI',32,9);
CPDI_displacement9(3) = results(1);CPDI_velocity9(3) = results(2);CPDI_stress9(3) = results(3);
clear results
results = Vibration('CPDI',64,9);
CPDI_displacement9(4) = results(1);CPDI_velocity9(4) = results(2);CPDI_stress9(4) = results(3);
clear results
results = Vibration('CPDI',128,9);
CPDI_displacement9(5) = results(1);CPDI_velocity9(5) = results(2);CPDI_stress9(5) = results(3);
clear results
results = Vibration('CPDI',256,9);
CPDI_displacement9(6) = results(1);CPDI_velocity9(6) = results(2);CPDI_stress9(6) = results(3);
clear results
results = Vibration('CPDI',512,9);
CPDI_displacement9(7) = results(1);CPDI_velocity9(7) = results(2);CPDI_stress9(7) = results(3);
clear results
% % results = Vibration('CPDI',512*2);
% % CPDI_displacement(8) = results(1);CPDI_velocity(8) = results(2);CPDI_stress(8) = results(3);
% % clear results
% % results = Vibration('CPDI',512*4);
% % CPDI_displacement(9) = results(1);CPDI_velocity(9) = results(2);CPDI_stress(9) = results(3);
% % clear results
% % results = Vibration('CPDI',512*8);
% % CPDI_displacement(10) = results(1);CPDI_velocity(10) = results(2);CPDI_stress(10) = results(3);
% % clear results
% % 
%% IMLS_CPDI interpolator
% results = Vibration('IMLS_CPDI',8,1);
% IMLS_CPDI_displacement(1) = results(1);IMLS_CPDI_velocity(1) = results(2);IMLS_CPDI_stress(1) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',16,1);
% IMLS_CPDI_displacement(2) = results(1);IMLS_CPDI_velocity(2) = results(2);IMLS_CPDI_stress(2) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',32,1);
% IMLS_CPDI_displacement(3) = results(1);IMLS_CPDI_velocity(3) = results(2);IMLS_CPDI_stress(3) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',64,1);
% IMLS_CPDI_displacement(4) = results(1);IMLS_CPDI_velocity(4) = results(2);IMLS_CPDI_stress(4) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',128,1);
% IMLS_CPDI_displacement(5) = results(1);IMLS_CPDI_velocity(5) = results(2);IMLS_CPDI_stress(5) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',256,1);
% IMLS_CPDI_displacement(6) = results(1);IMLS_CPDI_velocity(6) = results(2);IMLS_CPDI_stress(6) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',512,1);
% IMLS_CPDI_displacement(7) = results(1);IMLS_CPDI_velocity(7) = results(2);IMLS_CPDI_stress(7) = results(3);
% clear results
% 
% results = Vibration('IMLS_CPDI',8,4);
% IMLS_CPDI_displacement4(1) = results(1);IMLS_CPDI_velocity4(1) = results(2);IMLS_CPDI_stress4(1) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',16,4);
% IMLS_CPDI_displacement4(2) = results(1);IMLS_CPDI_velocity4(2) = results(2);IMLS_CPDI_stress4(2) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',32,4);
% IMLS_CPDI_displacement4(3) = results(1);IMLS_CPDI_velocity4(3) = results(2);IMLS_CPDI_stress4(3) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',64,4);
% IMLS_CPDI_displacement4(4) = results(1);IMLS_CPDI_velocity4(4) = results(2);IMLS_CPDI_stress4(4) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',128,4);
% IMLS_CPDI_displacement4(5) = results(1);IMLS_CPDI_velocity4(5) = results(2);IMLS_CPDI_stress4(5) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',256,4);
% IMLS_CPDI_displacement4(6) = results(1);IMLS_CPDI_velocity4(6) = results(2);IMLS_CPDI_stress4(6) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',512,4);
% IMLS_CPDI_displacement4(7) = results(1);IMLS_CPDI_velocity4(7) = results(2);IMLS_CPDI_stress4(7) = results(3);
% clear results
% 
% results = Vibration('IMLS_CPDI',8,9);
% IMLS_CPDI_displacement9(1) = results(1);IMLS_CPDI_velocity9(1) = results(2);IMLS_CPDI_stress9(1) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',16,9);
% IMLS_CPDI_displacement9(2) = results(1);IMLS_CPDI_velocity9(2) = results(2);IMLS_CPDI_stress9(2) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',32,9);
% IMLS_CPDI_displacement9(3) = results(1);IMLS_CPDI_velocity9(3) = results(2);IMLS_CPDI_stress9(3) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',64,9);
% IMLS_CPDI_displacement9(4) = results(1);IMLS_CPDI_velocity9(4) = results(2);IMLS_CPDI_stress9(4) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',128,9);
% IMLS_CPDI_displacement9(5) = results(1);IMLS_CPDI_velocity9(5) = results(2);IMLS_CPDI_stress9(5) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',256,9);
% IMLS_CPDI_displacement9(6) = results(1);IMLS_CPDI_velocity9(6) = results(2);IMLS_CPDI_stress9(6) = results(3);
% clear results
% results = Vibration('IMLS_CPDI',512,9);
% IMLS_CPDI_displacement9(7) = results(1);IMLS_CPDI_velocity9(7) = results(2);IMLS_CPDI_stress9(7) = results(3);
% clear results
% % % results = Vibration('IMLS_CPDI',512*2);
% % % IMLS_CPDI_displacement(8) = results(1);IMLS_CPDI_velocity(8) = results(2);IMLS_CPDI_stress(8) = results(3);
% % % clear results
% % % results = Vibration('IMLS_CPDI',512*4);
% % % IMLS_CPDI_displacement(9) = results(1);IMLS_CPDI_velocity(9) = results(2);IMLS_CPDI_stress(9) = results(3);
% % % clear results
% % % results = Vibration('IMLS_CPDI',512*8);
% % % IMLS_CPDI_displacement(10) = results(1);IMLS_CPDI_velocity(10) = results(2);IMLS_CPDI_stress(10) = results(3);
% % % clear results

%% Plot results
% dcell = [0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625 0.001953125 0.001953125/2 0.001953125/4 0.001953125/8];
dcell = [0.125 0.0625 0.03125 0.015625 0.0078125 0.00390625 0.001953125];
% dcell = [1 5 10 20 50 100];

figure 
loglog(dcell,MPM_displacement,':ro',dcell,CPDI_displacement,':dg',dcell,CPDI_displacement4,'--dg',dcell,CPDI_displacement9,'-dg',dcell,IMLS_CPDI_displacement,'-sb','linewidth',3)
legend('MPM-1PPC','CPDI-1PPC','CPDI-4PPC','CPDI-9PPC','CPLS-1PPC','Location','northeast')
xlabel('cell size (m)')
ylabel('Time (s)')
set(gca,'fontsize', 15)


figure 
% loglog(dcell,MPM_displacement,':ro',dcell,CPDI_displacement,'--dg',dcell,IMLS_CPDI_displacement,'-sb','linewidth',3)
loglog(dcell,CPDI_displacement,':dg',dcell,IMLS_CPDI_displacement,':sb','linewidth',3)
hold on
loglog(dcell,CPDI_displacement4,'--dg',dcell,IMLS_CPDI_displacement4,'--sb','linewidth',3)
hold on
loglog(dcell,MPM_displacement9,'-ro',dcell,CPDI_displacement9,'-dg',dcell,IMLS_CPDI_displacement9,'-sb','linewidth',3)
hold on
x1 = [0.002 0.2]; y1 = [0.00002 0.2];
% x1 = [0.002 0.2]; y1 = [0.00000002 0.0002];
plot(x1,y1,'--k','linewidth',2);
hold on
% legend('MPM','CPDI','CPLS','2nd order','Location','southeast')
string = {'CPDI-1PPC';'CPLS-1PPC';'CPDI-4PPC';'CPLS-4PPC';'MPM-9PPC';'CPDI-9PPC';'CPLS-9PPC';'2nd order'};
columnlegend(2,string,'Location','southeast','fontsize', 15)
xlabel('cell size (m)')
ylabel('RMS')
set(gca,'fontsize', 15)
axis([-inf inf 0.0000001 0.1]);

figure 
% loglog(dcell,MPM_velocity,':ro',dcell,CPDI_velocity,'--dg',dcell,IMLS_CPDI_velocity,'-sb','linewidth',3)
loglog(dcell,CPDI_velocity,':dg',dcell,IMLS_CPDI_velocity,':sb','linewidth',3)
hold on
loglog(dcell,CPDI_velocity4,'--dg',dcell,IMLS_CPDI_velocity4,'--sb','linewidth',3)
hold on
loglog(dcell,MPM_velocity9,'-ro',dcell,CPDI_velocity9,'-dg',dcell,IMLS_CPDI_velocity9,'-sb','linewidth',3)
hold on
% x2 = [0.002 0.2]; y2 = [0.0001 1];
x2 = [0.002 0.2]; y2 = [0.01 100];
plot(x2,y2,'--k','linewidth',2);
hold on
% legend('MPM','CPDI','CPLS','2nd order','Location','southeast')
string = {'CPDI-1PPC';'CPLS-1PPC';'CPDI-4PPC';'CPLS-4PPC';'MPM-9PPC';'CPDI-9PPC';'CPLS-9PPC';'2nd order'};
columnlegend(2,string,'Location','southeast','fontsize', 15)
xlabel('cell size (m)')
ylabel('RMS')
set(gca,'fontsize', 15)
axis([-inf inf 0.00001 inf]);


figure 
% loglog(dcell,MPM_stress,':ro',dcell,CPDI_stress,'--dg',dcell,IMLS_CPDI_stress,'-sb','linewidth',3)
loglog(dcell,CPDI_stress,':dg',dcell,IMLS_CPDI_stress,':sb','linewidth',3)
hold on
loglog(dcell,CPDI_stress4,'--dg',dcell,IMLS_CPDI_stress4,'--sb','linewidth',3)
hold on
loglog(dcell,MPM_stress9,'-ro',dcell,CPDI_stress9,'-dg',dcell,IMLS_CPDI_stress9,'-sb','linewidth',3)
hold on
x3 = [0.002 0.2]; y3 = [10 100000];
plot(x3,y3,'--k','linewidth',2);
hold on
% legend('MPM','CPDI','CPLS','2nd order','Location','southeast')
string = {'CPDI-1PPC';'CPLS-1PPC';'CPDI-4PPC';'CPLS-4PPC';'MPM-9PPC';'CPDI-9PPC';'CPLS-9PPC';'2nd order'};
columnlegend(2,string,'Location','southeast','fontsize', 15)
xlabel('cell size (m)')
ylabel('RMS')
set(gca,'fontsize', 15)
axis([-inf inf 1 inf]);
