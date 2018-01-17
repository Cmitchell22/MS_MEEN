%% Startup Algorithm for FDE
clear all
close all
clc
%% Initial Values
alpha = 0.25;
fdefun = @(t,y) (40320/gamma(9-alpha))*t^(8-alpha) - 3*(gamma(5+alpha/2)/gamma(5-alpha/2))*t^(4-alpha/2) + 9/4*gamma(alpha+1) + (3/2*t^(alpha/2)-t^4)^3-y^(3/2);
h = 0.05;
g = zeros(1,10)
f_corrector(1) = feval(fdefun, 0, 0);

%% Approximate y(h/4) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predict y(h/4):
jj = 1
y_predictor(jj+1) = g(jj+1) + 1/gamma(alpha + 1)*(h/4)^alpha*f_corrector(jj)

%% Correct y(h/4)
%Intermediate Steps:
% 1. Calculate B coefficients
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(0, h/4, h/4, alpha) 
% 2. Calculate f_pred, using y_pred.
f_predictor(jj+1) = feval(fdefun, h/4, y_predictor(jj+1)) 

y_corrector(jj+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*f_predictor(jj+1))

% Calculate f_corrector(h/4)
f_corrector(jj+1) = feval(fdefun, h/4, y_corrector(jj+1))

jj = jj + 1;
% jj + 1 = 3, corresponding to h/2.
%     jj = 2, corresponding to h/4.
% jj - 1 = 1, corresponding to 0.
%% Approximate y(h/2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lag term:
%Intermediate Steps:
%1. Calculate B Coefficients
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(0, h/4, h/2, alpha) 


y_lag(jj+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj-1) + B_1n_Coeff*f_corrector(jj))
%% Predict y(h/2) Constant Interp:
y_predictor_P1(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha+1)*(h/4)^alpha*f_corrector(jj)

%% Predict y(h/2) Linear Interp:
%Intermediate Steps:
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(h/4, h/2, h/2, alpha) 

y_predictor_P2(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha+1)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*y_predictor_P1(jj+1));

%% Correct y(h/2) : Quadratic Interp:
%Intermediate Steps:
[A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(0,h/2,0,h/4,h/2,h/2,alpha)

y_corrector(jj+1) = g(jj+1) + 1/gamma(alpha+1)*(A_0n_Coeff*f_corrector(jj-1) + A_1n_Coeff*y_corrector(jj) + A_2n_Coeff*y_predictor_P2(jj+1));

f_corrector(jj+1) = feval(fdefun, h/2, y_corrector(jj+1))

jj = jj +1;
% jj + 1 = 4, corresponding to h.
%     jj = 3, corresponding to h/2.
% jj - 1 = 2, corresponding to h/4.
% jj - 2 = 1, corresponding to 0.
%% Approximate y(h) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lag term:
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(0, h/2, h, alpha) 


y_lag(jj+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj-2) + B_1n_Coeff*f_corrector(jj))

%% Predict y(h) Constant Interp:
y_predictor_P1(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha+1)*(h/2)^alpha*f_corrector(jj)

%% Predict y(h) Linear Interp:
%Intermediate Steps:
[B_0n_Coeff, B_1n_Coeff] = B_Coefficients(h/2, h, h, alpha) 

y_predictor_P2(jj+1) = g(jj+1) + y_lag(jj+1) + 1/gamma(alpha+1)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*y_predictor_P1(jj+1));

%% Correct y(h) : Quadratic Interp:
%Intermediate Steps:
[A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(0,h,0,h/2,h,h,alpha)

y_corrector(jj+1) = g(jj+1) + 1/gamma(alpha+1)*(A_0n_Coeff*f_corrector(jj-2) + A_1n_Coeff*y_corrector(jj) + A_2n_Coeff*y_predictor_P2(jj+1));

Final_Y_Startup = y_corrector(jj+1)
jj = jj +1;
% jj + 1 = 4, corresponding to h.
%     jj = 3, corresponding to h/2.
% jj - 1 = 2, corresponding to h/4.
% jj - 2 = 1, corresponding to 0.

%% PROGRESS: SUCCESSFULLY CALCULATING Y(1) FOR ALGORITHM FROM STARTUP. NEED TO INTRODUCE INTITIAL Y_LAG AS WELL POSSIBLY?
