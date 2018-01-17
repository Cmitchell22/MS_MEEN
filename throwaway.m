y = @(x) x.^2
z = @(x) (x + 2)
A = @(x) y(x).*z(x)
quad(A,0,1)