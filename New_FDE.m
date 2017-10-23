syms t;
fdefun = t^2;
alpha  = 0.5;
b0 = -1;
b1 = alpha + 2;
y_lag_sum = 0; 
y_corrector_j = fdefun;
m = 1;
time_n = [0:0.1:1];
first_pass = true;
h = 0.1;


for n = 1:11
    
      for j = 1:n+1-1
      
        fprintf('making a pass\n');
        
        [B_0n_Coeff, B_1n_Coeff] = B_Coefficients((j-1)*h,(j+1-1)*h,(n+1-1)*h, alpha)
        
        f_corrector_j = Take_Caputo_Deriv(y_corrector_j, alpha, m)
        
        f_corrector_jj = Take_Caputo_Deriv(y_corrector_jj, alpha, m)
        
        y_lag_sum = y_lag_sum + B_0n_Coeff*f_corrector_j + B_1n_Coeff*f_corrector_jj
      end 
      
    g = calculate_IC(fdefun, alpha, n+1)
    
    y_lag(n+1) = 1/gamma(alpha)*y_lag_sum
    
    % Predictor term at the new step
    if n > 1
      y_predictor(n+1) = g + y_lag(n+1) + h^alpha/gamma(alpha+2)*(b0*f_corrector(n-1) + b1*f_corrector(n))
    else
        error('n is an indexing term; it cannot be <= 0, as Matlab is a 1-indexed program.')
    end
    
    f_predictor(n+1) = Take_Caputo_Deriv(y_predictor(n+1), alpha, m)
    
    f_corrector(n) = Take_Caputo_Deriv(y_corrector(n), alpha, m)
    
    if first_pass == false
        f_predictor(n+1) = Take_Caputo_Deriv(y_predictor(n+1), alpha, m)
    end
    
    [B_0n_Coeff, B_1n_Coeff] = B_Coefficients((n-1)*h,(n+1-1)*h, (n+1-1)*h, alpha)
    
    %Increment term at the new step
    y_increment(n+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(n) + B_1n_Coeff*f_predictor(n+1))

    %Corrector term at the new step
    y_corrector(n + 1) = g + y_lag(n + 1) + y_increment(n + 1)
     
    first_pass = false;
end