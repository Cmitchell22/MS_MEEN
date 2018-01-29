function [T,Y] = MF_PCQ(fdefun, t0, tfinal, y0, alpha, h)

%% Initialization of Variables, and I.C.'s
% alpha  = 1.25;
% y0 = [0,0];
% h = 0.05;
% t0 = 0;
% tfinal = 1;
%% Internal Variable Initialization
a0 = (alpha + 4)/2;
a1 = -2*(alpha+3);
a2 = (2*alpha^2 + 9*alpha + 12)/2;
y_lag_sum = 0;
m = 1; % TODO: should be calculated later with IC

%% Initialize Counters:
% NOTE: Since the step can be a non-integer value, (i.e. n =
% t0:h:tfinal), we need an integer variable corresponding to this step for
% array indexing purposes. These are defined below:
% 1) jj: This is the indexing variable for lag calculations. It corresponds
%        with the step j, where j is summed from 0 to n-1
% 2) nn: This is the indexing variable for the method. It corresponds with
% the step n, which is the current timestep. 
jj = 1;
nn = 3;

% Startup Export contains the values needed from the startupProcedure. As
% defined below, StartupProcedure contains:
% 1) startup_export(1) = f_corrector evaluated at t = h/2
% 2) startup_export(2) = f_corrector evaluated at t = h
% 3) startup_export(3) = f_corrector evaluated at t = 2h
% 4) startup_export(4) = y_corrector evaluated at t = h/2
% 5) startup_export(5) = y_corrector evaluated at t = h
% 6) startup_export(5) = y_corrector evaluated at t = 2h
scheme = 'quadratic';
startup_export = StartupProcedure(fdefun, scheme, t0, y0, alpha, h);
f_corrector(1) = feval(fdefun, t0, y0(1));
f_corrector(2) = startup_export(2);
f_corrector_one_half = startup_export(1);
f_corrector(3) = startup_export(3)

y_corrector(1) = y0(1);
y_corrector(2) = startup_export(5);
y_corrector_one_half = startup_export(4);
y_corrector(3) = startup_export(6)

%Calculate Initial Condition TODO *********
g_equation = 0; %calculate_IC(y_term, alpha)
%%
% Indexing is set up as follows:
% [ t = 0, t = 0.05 , t = 0.1, t = 0.15 ...]
% For n = 0, 1st pass:
% [ t = n, t = n+h ...]
% For n = 0.05=h, 1st pass:
% [ t = n-h, t = n, t = n+h]

%% 2nd Order, Linear Interp

% n is the current time step. Starts at 2*h, and increases by h each time.
% This is because the startup procedure calculates up through 2*h, so the
% next y_(n+1) term desired is y_(2*h + h). % This will calculate values up
% to final time tfinal, since each loop solves the y_(n+1) term.
for n = (2*h + t0):h:tfinal-h
    
      jj = 2;
      y_lag_sum = 0;
      
      for j = h:h:n-h
        %Calculate Coefficients for lag terms:
        [A_0j_Coeff, A_1j_Coeff, A_2j_Coeff] = A_Coefficients(j, j+h,j-h, j, j+h, n+h, alpha);
                
        % Summation Portion of the lag term. 
        y_lag_sum = y_lag_sum + A_0j_Coeff*f_corrector(jj-1) + A_1j_Coeff*f_corrector(jj) + A_2j_Coeff*f_corrector(jj+1);
        
        %Increment the interior counter
        jj = jj + 1;
      end 
      
    %Calculate the initial condition term **Still unsure about this term**  
    g(nn+1) = 0; % subs(g_equation, t, n+h);
    
    [A_00_Coeff, A_10_Coeff, A_20_Coeff] = A_Coefficients(0,h,0,h/2,h,n+h, alpha);
    
    A_initial = A_00_Coeff*f_corrector(1) + A_10_Coeff*f_corrector_one_half + A_20_Coeff*f_corrector(2); %****TODO
    
    % Lag term for this time step is now complete!
    y_lag_final(nn+1) = 1/gamma(alpha)*(A_initial + y_lag_sum);
    %Now we are able to calculate the predictor:
    
    % Predictor term at the new step (Need to have n-1 term. First pass: n = 0, n+1 = n+h, n-1 is not bound.)
    y_predictor(nn+1) = g(nn+1) + y_lag_final(nn+1) + h^alpha/gamma(alpha+3)*(a0*f_corrector(nn-2) + a1*f_corrector(nn-1) + a2*f_corrector(nn));
    
    %Solve the FDE function using this y_predicted value, at the current time
    %step.
    f_predictor(nn+1) = feval(fdefun, n+h, y_predictor(nn+1));
    
    %Calculate Coefficients for Linear Interp method, if not already
    %calculated. 
    if exist('A_0n_Coeff') == 0 || exist('A_1n_Coeff') == 0 || exist('A_2n_Coeff') == 0
        [A_0n_Coeff, A_1n_Coeff, A_2n_Coeff] = A_Coefficients(n, n+h, n-h, n, n+h, n+h, alpha);
    end
    
    %Increment term at the new step
    y_increment(nn+1) = 1/gamma(alpha)*(A_0n_Coeff*f_corrector(nn-1) + A_1n_Coeff*f_corrector(nn) + A_2n_Coeff*f_predictor(nn+1));

    %Corrector term at the new step
    y_corrector(nn+1) = g(nn+1) + y_lag_final(nn+1) + y_increment(nn+1);
      
    %Solve using corrector y_value guess instead of predicted y_value. To
    %be used in the next prediction stage.
    f_corrector(nn+1) = feval(fdefun, n+h, y_corrector(nn+1));
    
    nn = nn + 1;
     
end 
%returns
Y = y_corrector;
T = t0:h:tfinal;
%% END of Method
end
