function [L_j, L_jj] = Linear_Interp(t_j, t_jj)
%t_a = lower bound of integration
%t_b = upper bound of integration
%t_c = 
syms tau
L_j  = (tau - t_jj) / (t_j - t_jj);

L_jj = (tau - t_j) / (t_jj - t_j);

end
