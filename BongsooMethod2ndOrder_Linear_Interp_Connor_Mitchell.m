function [t, y] = fde_BongsooJang(alpha, fdefun, t0, tfinal, y0, h)
%Insert Help Description Here

%% Error Handling %%
if nargin < 6
	error('This Function requires 6 inputs: alpha, fdefun, t0, tfinal, y0, and h')
end

if alpha <= 0
    error('The order alpha or the FDE must be positive.') 
end

%Expand Error Handling!!!

%% variable initialization and initial values

SUMMER_1 = 0;
b_0n = -1;
b_1n = alpha + 2;
n = (tfinal-t0)/h;
y(0) = y0;
STEP = 0.25;
g(0) = 0;

%Logic for initialization. the first prediction will be y(1/4), then y(1/2), then y(1), and all further steps are 1.
 
for t_j = t0:STEP:n-1  %timestep. t_j is the current time t(j), t_jj is t(j+1), the next step.

	if t_j == t0

	    STEP = 0.25;

	elseif t_j == t0 + 0.25

		STEP = 0.5;

	else STEP = 1;
        
 t_jj = t_j + STEP;

	end

     %%initial condition term
	 g(t_jj) = g(t_j) + y(t_j)*(t_jj)^(t_j)/factorial(t_j) %***************

    % Constant Interpolation
    

    %This L function is the linear interpolator, as definded in Equation (2.2)
    L(t_j)  = @(tau) (tau - t_jj) / (t_j - t_jj);
    L(t_jj) = @(tau) (tau - t_j) / (t_jj - t_j);

    %Under Lemma 2.1 There are two different forms for calculating B coefficients, both are presented below. 
    %The latter is a more specific case of use.
    % This B_0n_Coeff and B_1n_Coeff are written in the form as expressed by Equations (2.5, 2.6),
    % Where phi(tau) is the Linear Interpolator in this case.
    if t_jj <= 1

    	B_0n_funct = @(tau) (t_jj-tau)^(alpha-1) * L(t_j)
    	B_0n_Coeff = integral(B_0n_funct,t_j,t_jj)

    	B_1n_funct = @(tau) (t_jj-tau)^(alpha-1) * L(t_jj)
    	B_1n_Coeff = integral(B_1n_funct,t_j,t_jj)

    else 

    	B_0n_Coeff(t_jj) = h^alpha/(alpha*(alpha+1))*(b_0n*L(t_j-1)+b_1n*L(t_j));
    	B_1n_Coeff(t_jj) = h^alpha/(alpha*(alpha+1))*(b_0n*L(t_j)+b_1n*L(t_jj));

	end

    %%Summation Term for computing the Y_lag term
    SUMMER_1 = SUMMER_1 + B_0n_Coeff(t_jj)*f_corrector(t_j) + B_1n_Coeff(t_jj)*f_corrector(t_jj);

    %y_lag term computed.
    y_lag(t_jj) = 1/gamma(alpha)*SUMMER_1;
    
    f_corrector(t_j)  = f_predictor(t_j);
    f_predictor(t_jj) = Take_Caputo_Deriv(y, t, alpha, t_j, t_jj, m)

    %% Y_increment term computation)
    Y_increment(t_jj) = 1/gamma(alpha)*(B_0n_Coeff(t_jj)*f_corrector(t_j-1) + B_1n_Coeff*f_predictor(t_j));

    % y_corrector term computed using sum of g(n+1), y_lag(n+1), and Y_increment(n+1), Equation (2.9)
    y_corrector(t_jj) = g(t_jj) + y_lag(t_jj) + Y_increment(t_jj);

    %Now the y_prediction(n+1)
    y_prediction(t_jj) = g(t_jj) + y_lag(t_jj) + h^alpha/gamma(alpha + 2)*(b_0n*f_corrector(t_j-1)+b_1n*f_corrector(t_j));
end