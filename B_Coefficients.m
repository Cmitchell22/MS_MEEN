function [B_0n_Coeff, B_1n_Coeff] = B_Coefficients(t_j, t_jj, t_nn, alpha)
%t_a = lower bound of integration
%t_b = upper bound of integration
syms tau
[L_j, L_jj] = Linear_Interp(t_j,t_jj);
B_0n_funct = ((t_nn-tau).^(alpha-1))*L_j;
 B_0n_Coeff = int(B_0n_funct,t_j,t_jj);
 
 B_1n_funct = ((t_nn-tau).^(alpha-1))*L_jj;
 B_1n_Coeff = int(B_1n_funct,t_j,t_jj);
end