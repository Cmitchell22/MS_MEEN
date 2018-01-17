function [y_deriv] = Take_Deriv(y,m)
if y == 0
    y_deriv = 0;
elseif m == 0
    y_deriv = y;
else
    y_deriv = diff(y,m);
end