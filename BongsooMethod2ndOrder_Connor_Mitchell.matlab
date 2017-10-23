%%Variable Initialization
t_j;
t_jj;
t_nn; % t_nn represents the t_n+1 time.
tau;
alpha;


B_nn0(0) = 0;
B_nn1(0) = 0;

for j = 0:n-1

	B_nn0(j+1) = integral( (t_nn-tau)^(alpha-1)*L_j(tau), j, j+1 );     % B_n+1, 0,j
    B_nn1(j+1) = integral( (t_nn-tau)^(alpha-1)*L_jj(tau), j, j+1 );    % B_n+1, 1,j
end

for j = 0:n-1
	SUMMER_1 = B_nn0(j)*f_corrector(j) + B_nn1*f_corrector(j+1)
end




y_corrector(t_nn) = g(t_nn) + 1/gamma(alpha)*SUMMER_1 ...
                  + 1/gamma(alpha) * (B_nn0_n*f_corrector(t_n)+B_nn1_n*f_predictor(t_nn);