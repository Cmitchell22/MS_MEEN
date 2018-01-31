function [startup_export] = Startup_Procedure(fdefun, scheme, t0, y0, alpha, h)

%% Initial Values
% alpha = 0.25;
% fdefun = @(t,y) (40320/gamma(9-alpha))*t^(8-alpha) - 3*(gamma(5+alpha/2)/gamma(5-alpha/2))*t^(4-alpha/2) + 9/4*gamma(alpha+1) + (3/2*t^(alpha/2)-t^4)^3-y^(3/2);
% h = 0.05;
g = zeros(1,10);
f_corrector(1) = feval(fdefun, 0, 0);

%% Approximate y(h/4) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predict y(h/4):
jj = 1;
y_predictor(jj+1) = g(jj+1) + 1/gamma(alpha + 1)*(h/4)^alpha*f_corrector(jj);

%% Correct y(h/4)
%Intermediate Steps:
% 1. Calculate B coefficients
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(0, h/4, h/4, alpha) ;
% 2. Calculate f_pred, using y_pred.
f_predictor(jj+1) = feval(fdefun, h/4, y_predictor(jj+1)) ;

y_corrector(jj+1) = g(jj+1) + 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*f_predictor(jj+1));

% Calculate f_corrector(h/4)
f_corrector(jj+1) = feval(fdefun, h/4, y_corrector(jj+1));

jj = jj + 1;
% jj + 1 = 3, corresponding to h/2.
%     jj = 2, corresponding to h/4.
% jj - 1 = 1, corresponding to 0.
%% Approximate y(h/2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lag term:
%Intermediate Steps:
%1. Calculate B Coefficients
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(0, h/4, h/2, alpha) ;


y_lag(jj+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj-1) + B_1n_Coeff*f_corrector(jj));
%% Predict y(h/2) Constant Interp:
y_predictor_P1(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha+1)*(h/4)^alpha*f_corrector(jj);

%% Predict y(h/2) Linear Interp:
%Intermediate Steps:
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(h/4, h/2, h/2, alpha) ;

f_predictor_P1(jj+1) = feval(fdefun, h/2, y_predictor_P1(jj+1));

y_predictor_P2(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*f_predictor_P1(jj+1));

%% Correct y(h/2) : Quadratic Interp:
%Intermediate Steps:
[A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(0,h/2,0,h/4,h/2,h/2,alpha);

f_predictor_P2(jj+1) = feval(fdefun, h/2, y_predictor_P2(jj+1));

y_corrector(jj+1) = g(jj+1) + 1/gamma(alpha)*(A_0n_Coeff*f_corrector(jj-1) + A_1n_Coeff*f_corrector(jj) + A_2n_Coeff*f_predictor_P2(jj+1));

f_corrector(jj+1) = feval(fdefun, h/2, y_corrector(jj+1));

jj = jj +1;
% jj + 1 = 4, corresponding to h.
%     jj = 3, corresponding to h/2.
% jj - 1 = 2, corresponding to h/4.
% jj - 2 = 1, corresponding to 0.
%% Approximate y(h) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lag term:
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(0, h/2, h, alpha) ;


y_lag(jj+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj-2) + B_1n_Coeff*f_corrector(jj));

%% Predict y(h) Constant Interp:
y_predictor_P1(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha+1)*(h/2)^alpha*f_corrector(jj);

%% Predict y(h) Linear Interp:
%Intermediate Steps:
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(h/2, h, h, alpha) ;

f_predictor_P1(jj+1) = feval(fdefun, h, y_predictor_P1(jj+1));


y_predictor_P2(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*f_predictor_P1(jj+1));

%% Correct y(h) : Quadratic Interp:
%Intermediate Steps:
[A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(0,h,0,h/2,h,h,alpha);

f_predictor_P2(jj+1) = feval(fdefun, h, y_predictor_P2(jj+1));


y_corrector(jj+1) = g(jj+1) + 1/gamma(alpha)*(A_0n_Coeff*f_corrector(jj-2) + A_1n_Coeff*f_corrector(jj) + A_2n_Coeff*f_predictor_P2(jj+1));

f_corrector(jj+1) = feval(fdefun, h, y_corrector(jj+1));

jj = jj +1;
% jj + 1 = 5, corresponding to 2h.
%     jj = 4, corresponding to h.
% jj - 1 = 3, corresponding to h/2.
% jj - 2 = 2, corresponding to h/4.
% jj - 3 = 1, corresponding to 0.
%% Approximate y(2h) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lag term:
[A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(0, h, 0, h/2, h, 2*h, alpha) ;

y_lag(jj+1) = 1/gamma(alpha)*(A_0n_Coeff*f_corrector(jj-3) + A_1n_Coeff*f_corrector(jj-1) + A_2n_Coeff*f_corrector(jj));

%% Predict y(2h) Constant Interp:
y_predictor_P1(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha+1)*(h)^alpha*f_corrector(jj);

%% Predict y(2h) Linear Interp:
%Intermediate Steps:
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(h, 2*h, 2*h, alpha) ;

f_predictor_P1(jj+1) = feval(fdefun, 2*h, y_predictor_P1(jj+1));


y_predictor_P2(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*f_predictor_P1(jj+1));

%% Correct y(2h) : Quadratic Interp:
%Intermediate Steps:
[A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(h,2*h,0,h,2*h,2*h,alpha);

f_predictor_P2(jj+1) = feval(fdefun, 2*h, y_predictor_P2(jj+1));


y_corrector(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha)*(A_0n_Coeff*f_corrector(jj-3) + A_1n_Coeff*f_corrector(jj) + A_2n_Coeff*f_predictor_P2(jj+1));

f_corrector(jj+1) = feval(fdefun, 2*h, y_corrector(jj+1));

% Linear vs Quadratic Startup:
if strcmp(scheme, 'linear') == true || strcmp(scheme, 'ABM') == true
    startup_export = [f_corrector(jj), f_corrector(jj+1), y_corrector(jj), y_corrector(jj+1)];
elseif strcmp(scheme, 'quadratic') == true
    startup_export = [f_corrector(jj-1), f_corrector(jj), f_corrector(jj+1), y_corrector(jj-1),  y_corrector(jj), y_corrector(jj+1)];
else
    error('Please enter a valid scheme for the startup algorithm to export the required values. The accepted inputs are "linear" and "quadratic".')
end
end

