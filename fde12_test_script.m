alpha = 0.25;
fdefun = @(t,y) (40320/gamma(9-alpha))*t^(8-alpha) - 3*(gamma(5+alpha/2)/gamma(5-alpha/2))*t^(4-alpha/2) + 9/4*gamma(alpha+1) + (3/2*t^(alpha/2)-t^4)^3-y^(3/2);
t0 = 0;
t_final = 1;
y0 = [0,0];
tfinal = 1;
h = 0.05;

[T,Y] = fde12(alpha,fdefun, t0, t_final,y0,h)