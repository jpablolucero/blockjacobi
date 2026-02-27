function [c0,u_exact,BC,rhs,poincareSteklovOperator] = ItITest1(k)

bump = @(x,y) -1.5 * exp(-160. * ((x-0.5).^2 + (y-0.5).^2));
c0 = @(x,y) -k^2*(1.0-bump(x,y));
u_exact = @(x, y) zeros(size(x));
du_exactdx = @(x, y) zeros(size(x));
du_exactdy = @(x, y) zeros(size(x));
BC = {@(x, y) - du_exactdx(x,y) + 1i * k .* u_exact(x,y),...
      @(x, y)   du_exactdx(x,y) + 1i * k .* u_exact(x,y),...
      @(x, y) - du_exactdy(x,y) + 1i * k .* u_exact(x,y),...
      @(x, y)   du_exactdy(x,y) + 1i * k .* u_exact(x,y)};
rhs = @(x, y) -k^2*bump(x,y).*exp(1i * k * x);

poincareSteklovOperator = "ItI";

end