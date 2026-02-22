function [u,A,rhs,px,py,A0,rhs0] = get_fem(div,u_fun,f_fun, ax, bx, ay, by, tol, mx, my)
if nargin < 10, my  = 2^div; end
if nargin < 9,  mx  = 2^div; end
if nargin < 8,  tol = 1e-14; end
if nargin < 7,  by  = 1; end
if nargin < 6,  ay  = 0; end
if nargin < 5,  bx  = 1; end
if nargin < 4,  ax  = 0; end

hx = (bx - ax) / mx;
hy = (by - ay) / my;

[xg,yg] = meshgrid(ax:hx:bx, ay:hy:by);

Nn   = (mx+1)*(my+1);
node = @(i,j) (j-1)*(mx+1)+i;

K = sparse(Nn,Nn);
F = zeros(Nn,1);

ni = 6;
[t, wt] = lgwt(ni, -1, 1);
gp = 0.5*(t+1);
w  = 0.5*wt;

for j = 1:my
    for i = 1:mx
        nd = [node(i,j), node(i+1,j), node(i,j+1), node(i+1,j+1)];
        x0 = ax + (i-1)*hx;
        y0 = ay + (j-1)*hy;

        Ke = zeros(4,4);
        Fe = zeros(4,1);

        for a = 1:ni
            X = gp(a); wx = w(a);
            for b = 1:ni
                Y = gp(b); wy = w(b);

                N   = [(1-X)*(1-Y); X*(1-Y); (1-X)*Y; X*Y];
                dNx = [  -(1-Y)/hx;(1-Y)/hx;   -Y/hx;Y/hx];
                dNy = [  -(1-X)/hy;   -X/hy;(1-X)/hy;X/hy];

                wgt = hx*hy*wx*wy;

                Ke = Ke + (dNx*dNx' + dNy*dNy') * wgt;
                Fe = Fe + N * f_fun(x0 + hx*X, y0 + hy*Y) * wgt;
            end
        end

        K(nd,nd) = K(nd,nd) + Ke;
        F(nd)    = F(nd) + Fe;
    end
end

vec = @(M) reshape(M.',[],1);
g   = vec(u_fun(xg,yg));

B = vec( (abs(xg-ax) < tol) | (abs(xg-bx) < tol) | ...
    (abs(yg-ay) < tol) | (abs(yg-by) < tol) );
I = ~B;

A   = K;
rhs = F;
rhs(I) = rhs(I) - K(I,B) * g(B);

A0   = A;
rhs0 = rhs;

A(:,B) = 0;  A(B,:) = 0;  A(B,B) = speye(nnz(B));
rhs(B) = g(B);

u  = A \ rhs;

px = vec(xg);
py = vec(yg);
end
