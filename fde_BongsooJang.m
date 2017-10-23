%function [t, y] = fde_BongsooJang(alpha, fdefun, t0, t_final, y0, h)
%Insert Help Description Here

clc; clear all; close all;
syms 't';
alpha = 0.5;
fdefun = t^2;
t0 = 1;
t_final = 100;
y0 = 0;
h = 1;
%% ********************** Error Handling **********************
% if nargin < 6
% 	error('This Function requires 6 inputs: alpha, fdefun, t0, t_final, y0, and h')
% elseif isnumeric(alpha) == 0 || isnumeric(t0) == 0 || isnumeric(t_final) == 0 || isnumeric(y0) == 0 || isnumeric(h) == 0
%     error('One of the inputs is not of the correct data type.')
% elseif alpha <= 0
%     error('The order alpha or the FDE must be positive.') 
% elseif t0 < 0
%     error('The initial time specified, t0, is less than zero. t0 must be greater than or equal to zero')
% elseif t_final < t0
%     error('The final time specified, t_final, is less than the initial time, t0. Please input a correct t_final.')
% end


%********Expand Error Handling***************
%% ************ Time/Iteration Setup **************************
startup_time = [0, 0.25, 0.5, 1, 2];
len = length(startup_time);
syms 't'

%% ****** Variable Initialization and Initial Conditions ******

SUMMER_1 = 0;
n = (t_final-t0)/h;
g = 0;
 % f_corrector = zeros(1,len);
 f_corrector(1) = 0*t;
% f_predictor = zeros(1,len);
% y_predictor = zeros(1,len);
b0_coeff = -1;
b1_coeff = alpha + 2;
m = ceil(alpha); %next integer larger than alpha
if m == alpha
    m = alpha + 1;
end

%% Initial Condition, g term 


     %%initial condition term
     for k = 0:m
        g = g + Take_Deriv(fdefun,k)*(t)^(k)/factorial(k) %***************
     end

%% ************************************************************

%Logic for initialization. the first prediction will be y(1/4), then y(1/2), then y(1), and all further steps are 1.
 
% for t_step = 1:(t_final-t0 + length(startup_time))  %timestep. t_j is the current time t(j), t_jj is t(j+1), the next step.

for t_step = 1:4 %(steps required for startup sequence
    fprintf('The Current step count is: %d', t_step)
    t_j  = startup_time(t_step);
    t_jj = startup_time(t_step + 1);    
     
%% ************************* Startup Algorithm *************************


%StartUp
if t_step == 1
    %Approximate y(1/4):
    
    %Compute Interpolation B coefficients
    t_a = startup_time(1); 
    t_b = startup_time(t_step +1);
    t_c = t_a;
    t_d = t_b;
    [B_0n_Coeff, B_1n_Coeff] = Linear_Interp(t_a,t_b,t_c,t_d,t_jj,alpha);
    
    %Predict: y_predictor(1/4): constant interpolation
    y_predictor(t_step + 1) = g + 1/(gamma(alpha+1))*(h/4)^alpha * f_corrector(t_step);
    f_predictor(t_step + 1) = Take_Caputo_Deriv( y_predictor(t_step + 1), alpha, m);
    %Correct: y_corrector(1/4): linear interpolation
    y_corrector(t_step + 1) = g + (1/gamma(alpha)) * ( B_0n_Coeff*f_corrector(t_step) + B_1n_Coeff*f_predictor(t_step + 1))
    y_corrector(t_step + 1) = subs(y_corrector(t_step + 1),t,t_step +1)
elseif t_step == 2
    % Approximate y(1/2):
    
    %Calculate Interpolation A & B coefficients for lag term
    t_a = startup_time(t_step-1);
    t_b = startup_time(t_step);
    t_c = t_a;
    t_d = t_b;
    [B_0n_Coeff, B_1n_Coeff] = Linear_Interp(t_a,t_b,t_c,t_d,t_jj,alpha);
    
    %Lag Term: y_lag(1/2):
    f_corrector(t_step) = Take_Caputo_Deriv(y_corrector(t_step), alpha, m)
    y_lag(t_step + 1) = 1/(gamma(alpha))* ( B_0n_Coeff*f_corrector(t_step -1) +B_1n_Coeff*f_corrector(t_step) )

   
    %Predict: y_predictor(1/2): const interp
    y_predictor_const(t_step + 1) = g + y_lag(t_step + 1) + 1/(gamma(alpha+1))*(h/4)^alpha * f_corrector(t_step)
    
      %Redefine Interpolation Values for Linear Interp Predictor
      t_a = startup_time(t_step);
      t_b = startup_time(t_step + 1);
      t_c = t_a;
      t_d = t_b;
      [B_0n_Coeff, B_1n_Coeff] = Linear_Interp(t_a,t_b,t_c,t_d,t_jj,alpha);
    
    %Predict: y_predictor(1/2): linear interpolation
    f_predictor_const(t_step + 1) = Take_Caputo_Deriv( y_predictor_const(t_step + 1), alpha, m)
    y_predictor_P1(t_step + 1) = g + y_lag(t_step + 1) + 1/gamma(alpha) * ( B_0n_Coeff*f_corrector(t_step) + B_1n_Coeff*f_predictor_const(t_step + 1))
    
      %Redefine Interpolation Values for Qudratic Interp Corrector
      t_a = startup_time(1);
      t_b = startup_time(t_step + 1);
      t_c = t_a;
      t_d = startup_time(t_step);
      t_e = t_b;
      [A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = Quadratic_Interp(t_a,t_b,t_c,t_d,t_e,t_jj,alpha);
    
    %Correct y_corrector(1/2): quadratic interpolation
    f_predictor_lin(t_step + 1) = Take_Caputo_Deriv( y_predictor_P1(t_step + 1), alpha, m)
    y_corrector(t_step +1) = g + 1/gamma(alpha) * ( A_0n_Coeff*f_corrector(t_step -1) +A_1n_Coeff*f_corrector(t_step) + A_2n_Coeff*f_predictor_lin(t_step +1))

elseif t_step == 3
    %Approximate y(1):
    
    %Lag Term: y_lag(1)
    
    %Calculate Interpolation B coefficients for lag term
    t_a = startup_time(1);
    t_b = startup_time(t_step);
    t_c = t_a;
    t_d = t_b;
    [B_0n_Coeff, B_1n_Coeff] = Linear_Interp(t_a,t_b,t_c,t_d,t_jj,alpha);
    
    f_corrector(t_step) = Take_Caputo_Deriv(y_corrector(t_step), alpha, m)

    y_lag(t_step + 1) = g + 1/gamma(alpha) * (B_0n_Coeff*f_corrector(t_step - 2) + B_1n_Coeff*f_corrector(t_step));
    %Predict: y_predictor(1) const interp
    y_predictor_P1(t_step + 1) = g + y_lag(t_step + 1) + 1/(gamma(alpha+1))*(h/2)^alpha * f_corrector(t_step)
    
    %Redefine Interpolation Values for Linear Interp Predictor
      t_a = startup_time(t_step);
      t_b = startup_time(t_step + 1);
      t_c = t_a;
      t_d = t_b;
      [B_0n_Coeff, B_1n_Coeff] = Linear_Interp(t_a,t_b,t_c,t_d,t_jj,alpha);
    
    %Predict: y_predictor(1) linear interp
    f_predictor_P1(t_step + 1) = Take_Caputo_Deriv(y_predictor_P1(t_step + 1), alpha, m)

    y_predictor_P2(t_step + 1) = g + y_lag(t_step + 1) + 1/gamma(alpha) * ( B_0n_Coeff*f_corrector(t_step) + B_1n_Coeff*f_predictor_P1(t_step + 1))
  
    %Redefine Interpolation Values for Qudratic Interp Corrector
      t_a = startup_time(1);
      t_b = startup_time(t_step + 1);
      t_c = t_a;
      t_d = startup_time(t_step);
      t_e = t_b;
      [A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = Quadratic_Interp(t_a,t_b,t_c,t_d,t_e,t_jj,alpha);
    
      
    f_predictor_P2(t_step + 1) = Take_Caputo_Deriv(y_predictor_P2(t_step + 1), alpha, m)

      
    %Correct: y_corrector(1) quadratic interp
    y_corrector(t_step + 1) = g + 1/gamma(alpha) * ( A_0n_Coeff*f_corrector(t_step -2) +A_1n_Coeff*f_corrector(t_step) + A_2n_Coeff*f_predictor_P2(t_step +1))
  
    %end Startup Procedure
elseif t_step == 4
    %Approximate y(t_step, for all  t_step > 3)
 
    %Calculate Interpolation A coefficients for lag term
    t_a = startup_time(1);
    t_b = startup_time(t_step);
    t_c = t_a;
    t_d = startup_time(t_step - 1);
    t_e = t_b;
    [A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = Quadratic_Interp(t_a,t_b,t_c,t_d,t_e,t_jj,alpha);
    
    f_corrector(t_step) = Take_Caputo_Deriv(y_corrector(t_step), alpha, m)
 
    %Lag Term: y_lag(t_step + 1)
    y_lag(t_step + 1) = 1/gamma(alpha) * ( A_0n_Coeff*f_corrector(1) + A_1n_Coeff*f_corrector(3) + A_2n_Coeff*f_corrector(4) );
    %Predict: y_predictor(1) const interp
    y_predictor_P1(t_step + 1) = g + y_lag(t_step + 1) + 1/(gamma(alpha+1))*h^alpha * f_corrector(t_step)
    
    %Redefine Interpolation Values for Linear Interp Predictor
      t_a = startup_time(t_step);
      t_b = startup_time(t_step + 1);
      t_c = t_a;
      t_d = t_b;
      [B_0n_Coeff, B_1n_Coeff] = Linear_Interp(t_a,t_b,t_c,t_d,t_jj,alpha);
    
      
    f_predictor_P1(t_step+1) = Take_Caputo_Deriv(y_predictor_P1(t_step+1), alpha, m)
  
    %Predict: y_predictor(1) linear interp
    y_predictor_P2(t_step + 1) = g + y_lag(t_step + 1) + 1/gamma(alpha) * ( B_0n_Coeff*f_corrector(t_step) + B_1n_Coeff*f_predictor_P1(t_step + 1))
  
    
    %Redefine Interpolation Values for Qudratic Interp Corrector
      t_a = startup_time(t_step);
      t_b = startup_time(t_step + 1);
      t_c = startup_time(1);
      t_d = t_a;
      t_e = t_b;
     [A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = Quadratic_Interp(t_a,t_b,t_c,t_d,t_e,t_jj,alpha);
    
    f_predictor_P2(t_step+1) = Take_Caputo_Deriv(y_predictor_P2(t_step+1), alpha, m)
 
    %Correct: y_corrector(1) quadratic interp
    y_corrector(t_step + 1) = g + 1/gamma(alpha) * ( A_0n_Coeff*f_corrector(1) +A_1n_Coeff*f_corrector(t_step) + A_2n_Coeff*f_predictor_P2(t_step +1))
end

end





%% %%%%%%%%%%% Linear Interpolator %%%%%%%%%%%%%%%%%%%%%%
interpolator = 'none';
if strcmp(interpolator,'linear') == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Lag term coefficients and calculation
    if t_step > 1 % The first interpolation for lag term will be a null sum.
    t_a = Total_Time(t_step - 1);
    t_b = Total_Time(t_step);
    t_c = t_a;
    t_d = t_b;
    [B_0n_Coeff, B_1n_Coeff] = Linear_Interp(t_a,t_b,t_c,t_d,t_jj,alpha);
    
    SUMMER_1 = SUMMER_1 + B_0n_Coeff*f_corrector(t_step - 1) + B_1n_Coeff*f_corrector(t_step);
    
    y_lag(t_step + 1) = 1/gamma(alpha) * SUMMER_1 ;
    
    end
    
    %increment term coefficients and calculation
    t_a = Total_Time(t_step);
    t_b = Total_Time(t_step + 1);
    t_c = t_a;
    t_d = t_b;
    [B_0n_Coeff, B_1n_Coeff] = Linear_Interp(t_a,t_b,t_c,t_d,t_jj,alpha);

    y_increment(t_step + 1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(t_step) + B_1n_Coeff*f_predictor(t_step + 1));
   
    %Calculate Y Prediction and Correction values
    y_predictor(t_step + 1) = g + y_lag(t_step + 1) + (h^alpha)/gamma(alpha+2)*(b0_coeff*f_corrector(t_step -1) + b1_coeff*f_corrector(t_step));
    
    y_corrector(t_step + 1) = g + y_lag(t_step + 1) + Y_increment(t_step + 1);
    
    
end

%% %%%%%%%%%%% Quadratic Interpolator %%%%%%%%%%%%%%%%%%%%%%
if strcmp(interpolator,'quadratic') == 1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Quadratic Interpolator method goes here
else 
end

