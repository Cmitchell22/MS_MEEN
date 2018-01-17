function g_equation = calculate_IC(funct, alpha)
syms t
m = 1;
g_sum = 0;

for k = 0:m
    k
    g_sum = g_sum + Take_Deriv(funct,k)*(t)^(k)/factorial(k) %***************
end

g_equation = g_sum;