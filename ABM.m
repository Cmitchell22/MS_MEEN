%ABM
clc;
clear all;
%% Initialization of Variables, and I.C.'s
% syms t;
alpha  = 0.25;
b0 = -1;
b1 = alpha + 2; 
% y_term = t^(4+alpha)
% fdefun = gamma(5+alpha)/24*t^4 + t^(8+2*alpha) - y_term^2
% y_term = t^8 - 3*t^(4+alpha/2)+9/4*t^alpha;
fdefun = @(t,y) (40320/gamma(9-alpha))*t^(8-alpha) - 3*(gamma(5+alpha/2)/gamma(5-alpha/2))*t^(4-alpha/2) + 9/4*gamma(alpha+1) + (3/2*t^(alpha/2)-t^4)^3-y^(3/2);
y_lag_sum = 0; 
m = 1;
time_n = [0:0.1:1];
first_pass = true;
h = 0.05;
t0 = 0;
y0 = 0;
f_corrector(1) = feval(fdefun,t0,y0);

ii = 1;
jj = 1;

%Calculate Initial Condition
g_equation = 0; %calculate_IC(y_term, alpha)
%%
% Indexing is set up as follows:
% [ t = 0, t = 0.05 , t = 0.1, t = 0.15 ...]
% For n = 0, 1st pass:
% [ t = n, t = n+h ...]
% For n = 0.05=h, 1st pass:
% [ t = n-h, t = n, t = n+h]


%% ABM Method
% n is the current time step. Starts at zero, and increases by h each time.
for n = 0:h:1    
      ii = 1;
      y_lag_sum = 0;
      for j = 0:h:n-h
        
        %Calculate Coefficients for lag terms: CONFIRMED GOOD
        [B_0j_Coeff, B_1j_Coeff] = B_Coefficients(j,j+h, n+h, alpha);
                
        % Summation Portion of the lag term. 
        y_lag_sum = y_lag_sum + B_0j_Coeff*f_corrector(ii) + B_1j_Coeff*f_corrector(ii+1);
        
        %Increment the interior counter
        ii = ii + 1;
      end 
      
    %Calculate the initial condition term **Still unsure about this term**  
    g(jj+1) = 0; %subs(g_equation, t, n+h);
    
    % Lag term for this time step is now complete!
    y_lag_final(jj+1) = 1/gamma(alpha)*y_lag_sum;
    
    %Now we are able to calculate the predictor:
    y_predictor_sum = 0;
    mm = 1;
    for j = 0:h:n
    
        C_Coeff = h^alpha/alpha*((n+1-j)^alpha - (n-j)^alpha);
        y_predictor_sum = y_predictor_sum + C_Coeff*f_corrector(mm);
        mm = mm+1;
    end
    
    % Predictor term at the new step (Need to have n-1 term. First pass: n = 0, n+1 = n+h, n-1 is not bound.)
    y_predictor(jj+1) = g(jj+1) + 1/gamma(alpha)*C_Coeff*y_predictor_sum;

    % take Caputo Derivative of the predictor
    f_predictor(jj+1) = feval(fdefun,n+h, y_predictor(jj+1));
    
    %Calculate Coefficients for
    if exist('B_0n_Coeff') == 0 || exist('B_1n_Coeff') == 0
        [B_0n_Coeff, B_1n_Coeff] = B_Coefficients(n,n+h, n+h, alpha);
    end
    %Increment term at the new step
    y_increment(jj+1) = 1/gamma(alpha)*(B_0n_Coeff*f_corrector(jj) + B_1n_Coeff*f_predictor(jj+1));

    %Corrector term at the new step
    y_corrector(jj+1) = g(jj+1) + y_lag_final(jj+1) + y_increment(jj+1);
    
    f_corrector(jj+1) = feval(fdefun, n+h, y_corrector(jj+1));
    
    jj = jj + 1;
     
end

i = 1;
for t = 0:0.05:1
    true_solution(i) = t^8 - 3*t^(4+alpha/2)+9/4*t^alpha;
    i = i +1;
end

