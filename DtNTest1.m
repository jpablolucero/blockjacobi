function [u_exact,rhs,poincareSteklovOperator] = DtNTest1(k)

u_exact = @(x,y) sin(k*pi*x-0.1).*sin(k*pi*y-0.2);
rhs = @(x,y) 2*(k*pi)^2*sin(k*pi*x-0.1).*sin(k*pi*y-0.2);

poincareSteklovOperator = "DtN";

end