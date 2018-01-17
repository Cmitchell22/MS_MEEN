function [A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(integrand_low,integrand_high,t_low,t_mid,t_up,t_eval,alpha)
 
 [Q_low, Q_mid, Q_up] = Quadratic_Interp(t_low, t_mid, t_up);
 A_kernel = @(tau) ((t_eval-tau).^(alpha-1));
 A_0n_funct = @(tau) A_kernel(tau).*Q_low(tau);
 A_0n_Coeff = integral(A_0n_funct,integrand_low,integrand_high);

 A_1n_funct = @(tau) A_kernel(tau).*Q_mid(tau);
 A_1n_Coeff = integral(A_1n_funct,integrand_low,integrand_high);
        
 A_2n_funct = @(tau) A_kernel(tau).*Q_up(tau);
 A_2n_Coeff = integral(A_2n_funct,integrand_low,integrand_high);