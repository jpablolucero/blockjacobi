function [u_exact,rhs,c0,BC,poincareSteklovOperator] = DtNTest1(k)

c0 = @(x,y) zeros(size(x));
u_exact = @(x,y) sin(k*pi*x-0.1).*sin(k*pi*y-0.2);
BC = {@(x, y) u_exact(x,y),...
      @(x, y) u_exact(x,y),...
      @(x, y) u_exact(x,y),...
      @(x, y) u_exact(x,y)};
rhs = @(x,y) 2*(k*pi)^2*sin(k*pi*x-0.1).*sin(k*pi*y-0.2);

poincareSteklovOperator = "DtN";

end