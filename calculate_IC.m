function g = calculate_IC(fdefun, alpha, t_eval)

m = round(alpha, 0);
g_sum = 0;

for k = 0:m
    g_sum = g_sum + Take_Deriv(fdefun,k)*(t_eval)^(k)/factorial(k); %***************
end

g = subs(g_sum, t_eval);
