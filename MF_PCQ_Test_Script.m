% MF-PCQ Test Script
clc;
clear all;
close all;
alpha = 1.25;
y0 = [0,0];
t0 = 0;
tfinal = 1;
h = 0.05;
fdefun = @(t,y) (40320/gamma(9-alpha))*t^(8-alpha) - 3*(gamma(5+alpha/2)/gamma(5-alpha/2))*t^(4-alpha/2) + 9/4*gamma(alpha+1) + (3/2*t^(alpha/2)-t^4)^3-y^(3/2);
i = 1;
for t = 0:h:1
        true_solution(i) = t^8 - 3*t^(4+alpha/2)+9/4*t^alpha;
        i = i+1;
    end
    time = 0:h:1;
    plot(time, true_solution);
[T, Y_lin] = MF_PCL(fdefun, t0, tfinal, y0, alpha, h);
hold on
plot(time, Y_lin)
[T,Y_quad] = MF_PCQ(fdefun, t0, tfinal, y0, alpha, h);
hold on
plot(time, Y_quad)
[T,Y_ABM] = fde12(alpha, fdefun, t0, tfinal, y0, h);
hold on
plot(time, Y_ABM)
legend('True', 'Linear', 'Quadratic', 'ABM')
savename = strcat('Eq1_alpha=',num2str(alpha),'_N=',num2str(1/h),'.png');
saveas(gcf, savename)