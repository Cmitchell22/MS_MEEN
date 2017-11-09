function [A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(integrand_low,integrand_high,t_low,t_mid,t_up,t_eval,alpha)
 
 syms tau
 [Q_low, Q_mid, Q_up] = Quadratic_Interp(t_low, t_mid, t_up);
 A_0n_funct = ((t_eval-tau)^(alpha-1)) * Q_low;
 A_0n_Coeff = int(A_0n_funct,integrand_low,integrand_high);

 A_1n_funct = ((t_eval-tau)^(alpha-1)) * Q_mid;
 A_1n_Coeff = int(A_1n_funct,integrand_low,integrand_high);
        
 A_2n_funct = ((t_eval-tau)^(alpha-1)) * Q_up;
 A_2n_Coeff = int(A_2n_funct,integrand_low,integrand_high);