function y = f_Bruss(t,x,param)
a = param(1) ; mu = param(2);
y(1,1) = a -(mu+1)*x(1)+x(1)^2*x(2) ;
y(2,1) = mu*x(1)-x(1)^2*x(2) ;