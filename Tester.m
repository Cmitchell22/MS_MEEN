%% Startup
clc
close all
clear all

%% Parameter Definitions
alpha = 0.8;
a = 1; 
mu = 4;
t0 = 0;
tfinal = 100;
y0 = [0.2 ; 0.03];
h = 2^(-6);

%% FDEFUN Definition
fdefun = @(t,y) [a-(mu+1)*y(1)+y(1)^2*y(2) ; mu*y(1)-y(1)^2*y(2) ] ; 

%% Output
[t, y_predictor] = fde_BongsooJang(alpha,fdefun, t0,tfinal,y0,h);

%% Plotting

figure(1)
plot(t,y_predictor(1,:),t,y_predictor(2,:));
xlabel('t') ; ylabel('y(t)');
legend('y_1(t)','y_2(t)');
title('FDE solved by the FDE12.m code');

