function [c0,u_exact,BC,rhs,poincareSteklovOperator] = ItITest2(k)

bump = @(x,y) -1.5 * exp(-160. * ((x-0.5).^2 + (y-0.5).^2));
c0 = @(x,y) -k^2 * (1.0 - bump(x,y));
u_exact = @(x, y) exp(1i * k * x) .* cos(pi * y) + 1.;
du_exactdx = @(x, y) 1i * k * exp(1i * k * x) .* cos(pi * y);
du_exactdy = @(x, y) -pi * exp(1i * k * x) .* sin(pi * y);
BC = {@(x, y) - du_exactdx(x,y) + 1i * k .* u_exact(x,y),...
      @(x, y)   du_exactdx(x,y) + 1i * k .* u_exact(x,y),...
      @(x, y) - du_exactdy(x,y) + 1i * k .* u_exact(x,y),...
      @(x, y)   du_exactdy(x,y) + 1i * k .* u_exact(x,y)};
rhs = @(x, y) ( k^2 + (pi)^2 ) .* exp(1i * k * x) .* cos(pi * y) + c0(x,y) .* u_exact(x,y);

poincareSteklovOperator = "ItI";

end