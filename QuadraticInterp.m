function [A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = QuadraticInterp(t_a,t_b,t_c,t_d,t_e,t_jj,alpha)
 
 
 A_0n_funct = @(tau) ((t_jj-tau)^(alpha-1) *  ((tau - t_d) / (t_c - t_d)) * ((tau - t_e)/(t_c - t_e)));
 A_0n_Coeff = integral(A_0n_funct,t_a,t_b);

 A_1n_funct = @(tau) ((t_jj-tau)^(alpha-1) *  ((tau - t_c) / (t_d - t_c)) * ((tau - t_e)/(t_d - t_e)));
 A_1n_Coeff = integral(A_1n_funct,t_a,t_b);
        
 A_2n_funct = @(tau) ((t_jj-tau)^(alpha-1) *  ((tau - t_c) / (t_e - t_c)) * ((tau - t_d)/(t_e - t_d)));
 A_2n_Coeff = integral(A_2n_funct,t_a,t_b);