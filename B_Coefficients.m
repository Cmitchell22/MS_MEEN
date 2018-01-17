function [B_0n_Coeff, B_1n_Coeff] = B_Coefficients(t_j, t_jj, t_nn, alpha)
%t_a = lower bound of integration
%t_b = upper bound of integration
[L_j, L_jj] = Linear_Interp(t_j,t_jj);
B_kernel = @(tau) ((t_nn-tau).^(alpha-1));
B_0n_funct = @(tau) B_kernel(tau).*L_j(tau);
B_0n_Coeff = integral(B_0n_funct, t_j,t_jj);

 B_1n_funct = @(tau) B_kernel(tau).*L_jj(tau);
 B_1n_Coeff = integral(B_1n_funct, t_j,t_jj);

end