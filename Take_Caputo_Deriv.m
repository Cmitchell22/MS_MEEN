function [Caputo_Derivative] = Take_Caputo_Deriv(y,alpha,m)
% This function takes a fractional order derivative using the Caputo
% Derivative method. The inputs are as follows:
% y: the function to be evaluated
% t: the time variable
% alpha: The order of the derivative, i.e. d^(alpha)/dx^(alpha)
% t_j: The lower limit of evaluation
% t_jj: The Upper limit of evaluation
% m: the nearest integer greater than alpha
format short;
verif = true;
syms 't';
syms 'tau';

if m < alpha
    error('alpha can not be greater than m')
elseif (m-1) > alpha
    error('alpha must be larger than m-1')
end

Coeff1 = 1/gamma(m-alpha);
fun =(tau-t)^(m-1-alpha);
y_deriv = Take_Deriv(y,m);
Convolution = fun*y_deriv;
Caputo_Der = Coeff1*int(Convolution,t,0,tau);
Caputo_Deriv = subs(Caputo_Der,tau,t);
Caputo_Derivative = subs(Caputo_Deriv,t);


%% Verif:
if verif == true
    exponent = 3;
    verif_ans = gamma(exponent+1)/(gamma(exponent+1-alpha))*t^(exponent-alpha);
end

Diff = verif_ans - Caputo_Derivative;
