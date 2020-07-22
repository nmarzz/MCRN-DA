close all; clear all;clc;
%% L96
a=load('RMSE_l96_POD20.mat');%RMSE for L96 with POD r=20
b=load('RMSE_l96_POD30.mat');%RMSE for L96 with POD r=30
c=load('RMSE_l96_POD20DMD20.mat');%RMSE for L96 with POD and DMD r=20
d=load('RMSE_l96.mat');%I-projection
Time=a.Time;
RMSEsave_a=a.RMSEsave;
RMSEsave_proj_a=a.RMSEsave_proj;
RMSEsave_b=b.RMSEsave;
RMSEsave_proj_b=b.RMSEsave_proj;
RMSEsave_c=c.RMSEsave;
RMSEsave_proj_c=c.RMSEsave_proj;
RMSEsave_d=d.RMSEsave;


figure(1)
plot(Time,RMSEsave_a, 'c-.', 'LineWidth', 1)
grid on
hold on;
plot(Time,RMSEsave_proj_a,'m:','LineWidth', 1)
plot(Time,RMSEsave_d,'b:','LineWidth', 1)
title('POD Projection(r=20)')
xlabel('Time')
ylabel('RMSE')
ylim([0 0.3])
legend('RMSE Original','RMSE Projected','Identity Projection','Location', 'Best')
hold off

figure(2)
plot(Time,RMSEsave_b, 'c.-', 'LineWidth', 1)
grid on
hold on;
plot(Time,RMSEsave_proj_b,'m:','LineWidth', 1)
plot(Time,RMSEsave_d,'b:','LineWidth', 1)
title('POD Projection(r=30)')
xlabel('Time')
ylabel('RMSE')
ylim([0 0.3])
legend('RMSE Original','RMSE Projected','Identity Projection','Location', 'Best')
hold off

figure(3)
plot(Time,RMSEsave_c, 'c-.', 'LineWidth', 1)
grid on
hold on;
plot(Time,RMSEsave_proj_c,'m:','LineWidth', 1)
plot(Time,RMSEsave_d,'b:','LineWidth', 1)
title('POD and DMD Projection(r=20)')
xlabel('Time')
ylabel('RMSE')
ylim([0 0.3])
legend('RMSE Original','RMSE Projected','Identity Projection','Location', 'Best')
hold off
%% SWE
E=load('RMSE_POD10_1day.mat');%SWE for 1 day POD r=10
F=load('RMSE_POD20_1day.mat');%SWE for 1 day POD r=20
G=load('RMSE_NOproj_1day.mat');%SWE for 1 day POD r=30
H=load('RMSE_POD30_1day.mat'); %I-proj
Time_SWE=E.Time;
RMSEsave_E=E.RMSEsave;
RMSEsave_proj_E=E.RMSEsave_proj;
RMSEsave_F=F.RMSEsave;
RMSEsave_proj_F=F.RMSEsave_proj;
RMSEsave_G=G.RMSEsave;
RMSEsave_H=H.RMSEsave;
RMSEsave_proj_H=H.RMSEsave_proj;
figure(4)
plot(Time_SWE,RMSEsave_E, 'g-.', 'LineWidth', 2)
grid on
hold on;
plot(Time_SWE,RMSEsave_proj_E,'r:','LineWidth', 2)
plot(Time_SWE,RMSEsave_G,'b-','LineWidth', 2)
title('POD Projection(r=10)')
xlabel('Time')
ylabel('RMSE')
ylim([0 12])
legend('RMSE Original','RMSE Projected','Identity Projection','Location', 'Best')
hold off

figure(5)
plot(Time_SWE,RMSEsave_F, 'g-.', 'LineWidth', 2)
grid on
hold on;
plot(Time_SWE,RMSEsave_proj_F,'r:','LineWidth', 2)
plot(Time_SWE,RMSEsave_G,'b-','LineWidth', 2)
title('POD Projection(r=20)')
xlabel('Time')
ylabel('RMSE')
ylim([0 12])
legend('RMSE Original','RMSE Projected','Identity Projection','Location', 'Best')
hold off

figure(6)
plot(Time_SWE,RMSEsave_H, 'g-.', 'LineWidth', 2)
grid on
hold on;
plot(Time_SWE,RMSEsave_proj_H,'r:','LineWidth', 2)
plot(Time_SWE,RMSEsave_G,'b-','LineWidth', 2)
title('POD Projection(r=30)')
xlabel('Time')
ylabel('RMSE')
ylim([0 12])
legend('RMSE Original','RMSE Projected','Identity Projection','Location', 'Best')
hold off