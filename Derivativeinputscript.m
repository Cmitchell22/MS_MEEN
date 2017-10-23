syms t;
y = 2*t^3 - 5*t^2 + 3;
m = 1;
t_j = 0;
t_jj = 1;
alpha = 0.5;
Take_Caputo_Deriv(y,t,alpha,t_j,t_jj,m);