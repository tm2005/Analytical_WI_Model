close all
clear all
clc
%% crystals geometry
R0 = 100;
L0 = 1;

%%
h = 10; % distance from the origin (or center of rotation)

angle_step=0.1; %rotation step from numerical rotation

[P_unique,rs_unique] = numerical_rotation(R0,L0,h,angle_step, 800);

rh0 = triangle_exact(rs_unique,h,L0,R0);
rh0_square = square_exact(rs_unique,h,L0,R0);
rh0_dirac = dirac_response_scaled(rs_unique,h,L0,R0);

P_unique = P_unique';


%% Plot

f1=figure('DefaultAxesFontSize',16);
set(gcf, 'renderer', 'painters','Position', [100 100 1200 900]);

plot(rs_unique,rh0,'-','LineWidth',1.6,'Color',[0,0,1]), hold on, 
plot(rs_unique,rh0_square,'-','LineWidth',1.6,'Color',[0.2 0.46 0.2]), plot(rs_unique,rh0_dirac,'LineWidth',1.6,'Color','#009999'), plot(rs_unique,P_unique,':','LineWidth',3,'Color',[1,0,0])
legend('Triangle approximation','Square approximation','Rotation of Dirac line','Numerical rotation','Location','northeast')
title('Approximative and numerical plot, R_0 = 100, L_0 = 1, h = 10', 'FontSize', 16);
xlabel('Distance (r)', 'FontSize', 14);
ylabel('Probability (P)', 'FontSize', 14);
xlim([0 R0]);
ylim([0 max([max(P_unique*6/5), max(rh0*6/5)])] )
xticks([0 5:5:100])
grid on

% saveas(f1,'Graf_svi_R0_100_L0_1_h_10.png')