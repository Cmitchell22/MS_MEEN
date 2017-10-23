function [y_deriv] = Take_Deriv(y,m)

syms t
y_deriv = diff(y,m)