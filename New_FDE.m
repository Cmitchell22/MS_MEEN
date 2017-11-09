clc;
clear all;
syms t;
fdefun = t^2;
alpha  = 0.5;
b0 = -1;
b1 = alpha + 2;
y_lag_sum = 0; 
m = 1;
time_n = [0:0.1:1];
first_pass = true;
h = 0.05;
y_corrector(1) = 0;

ii = 1;
jj = 1;

% Indexing is set up as follows:
% [ t = 0, t = 0.05 , t = 0.1, t = 0.15 ...]
% For n = 0, 1st pass:
% [ t = n, t = n+h ...]
% For n = 0.05=h, 1st pass:
% [ t = n-h, t = n, t = n+h]

% n is the current time step. Starts at zero, and increases by h each time.
for n = 0:h:h
    n
    
      ii = 1;
      y_lag_sum = 0;
      for j = 0:h:n-h
      j
        fprintf('making a pass\n');
        
        %Calculate Coefficients for lag terms: CONFIRMED GOOD
        [B_0j_Coeff, B_1j_Coeff] = B_Coefficients(j,j+h, n+h, alpha)
        
%         % Take Caputo Derivative to get f_tilde(t_j) from y_tilde(t_j)
%         f_corrector_j(ii) = Take_Caputo_Deriv(y_corrector(ii), alpha, m)
%         
%         % Take Caputo Derivative to get f_tilde(t_j+1) from y_tilde(t_j+1)
%         f_corrector_jj(ii+1) = Take_Caputo_Deriv(y_corrector(ii+1), alpha, m)
        
        % Summation Portion of the lag term. 
        y_lag_sum = y_lag_sum + B_0j_Coeff*f_corrector(ii) + B_1j_Coeff*f_corrector(ii+1)
        
        %Increment the interior counter
        ii = ii + 1;
      end 
      
    %Calculate the initial condition term **Still unsure about this term**  
    g(jj+1) = calculate_IC(fdefun, alpha, n+h)
    
    % Lag term for this time step is now complete!
    y_lag_final(jj+1) = 1/gamma(alpha)*y_lag_sum
    
    %Now we are able to calculate the predictor:
    
    % Take Caputo Derivatives of the y_corrector term at current step
    f_corrector(jj) = Take_Caputo_Deriv(y_corrector(jj),n, alpha, m)
    
    % Predictor term at the new step (Need to have n-1 term. First pass: n = 0, n+1 = n+h, n-1 is not bound.)
    if jj > 1
        y_predictor(jj+1) = g(jj+1) + y_lag_final(jj+1) + h^alpha/gamma(alpha+2)*(b0*f_corrector(jj-1) + b1*f_corrector(jj))
    else
        y_predictor(jj+1) = g(jj+1) + y_lag_final(jj+1) + h^alpha/gamma(alpha+2)*(b0*0 + b1*f_corrector(jj))
    end
    % take Caputo Derivative of the predictor
    f_predictor(jj+1) = Take_Caputo_Deriv(y_predictor(jj+1), n+h, alpha, m)
    
    %Calculate Coefficients for
    [B_0n_Coeff, B_1n_Coeff] = B_Coefficients(n,n+h, n+h, alpha)
    
    %Increment term at the new step
    y_increment(jj+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*f_predictor(jj+1))

    %Corrector term at the new step
    y_corrector(jj+1) = g(jj+1) + y_lag_final(jj+1) + y_increment(jj+1)
    
    f_corrector(jj+1) = Take_Caputo_Deriv(y_corrector(jj+1),n+h, alpha, m)
    
    jj = jj + 1
     
end