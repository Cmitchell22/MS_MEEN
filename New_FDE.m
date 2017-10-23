alpha  = 0.5;
b0 = -1;
b1 = alpha + 2;
y_lag_sum = 0; 
y_corrector(1) = 0
m = 1

n = 3;

    for j = 1:n-1
        
        [B_0n_Coeff, B_1n_Coeff] = B_Coefficients(j,j+1,n+1, alpha)
        
        f_corrector(j) = Take_Caputo_Deriv(y_corrector(j), alpha, m);
        f_corrector(j+1) = Take_Caputo_Deriv(y_corrector(j+1), alpha, m);
        
        y_lag_sum = y_lag_sum + B_0n_Coeff*f_corrector(j) + B_1n_Coeff*f_corrector(j+1)

    end 
    
    g = calculate_IC(fdefun, alpha, n+1)
    
    y_lag(n +1) = 1/gamma(alpha)*y_lag_sum;
    
    f_corrector(n) = Take_Caputo_Deriv(y_corrector(n), alpha, m);
    f_predictor(n+1) = Take_Caputo_Deriv(y_predictor(n+1), alpha, m);
    
    [B_0n_Coeff, B_1n_Coeff] = B_Coefficients(n,n+1,n+1, alpha)
    
    y_increment(n+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(n) + B_1n_Coeff*f_predictor(n+1));

    y_corrector(n + 1) = g + y_lag(n + 1) + y_increment(n + 1)

    y_predictor(n+1) = g + y_lag(n+1) + h^alpha/gamma(alpha+2)*(b0*f_corrector(n-1) + b1*f_corrector(n));
