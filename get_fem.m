function [px,py,K,rhs0,u] = get_fem(div,c0_fun,f_fun,ax,bx,ay,by,mx,my,BC)
    if nargin < 10
        BC = [];
    end
    if nargin < 9
        my = 2^div;
    end
    if nargin < 8
        mx = 2^div;
    end
    if nargin < 7
        by = 1;
    end
    if nargin < 6
        ay = 0;
    end
    if nargin < 5
        bx = 1;
    end
    if nargin < 4
        ax = 0;
    end

    u = [];

    hx = (bx - ax) / mx;
    hy = (by - ay) / my;

    [xg,yg] = meshgrid(ax:hx:bx, ay:hy:by);

    node = @(i,j) (j-1)*(mx+1)+i;
    Nn = (mx+1)*(my+1);

    K = sparse(Nn,Nn);
    rhs0 = zeros(Nn,1);

    quad = Quadrature(5);
    t = 0.5 * (quad.Points + 1);
    wt = 0.5 * quad.Weights;

    for j = 1:my
        for i = 1:mx
            nd = [node(i,j), node(i+1,j), node(i,j+1), node(i+1,j+1)];
            x0 = ax + (i-1)*hx;
            y0 = ay + (j-1)*hy;

            Ke = zeros(4,4);
            Fe = zeros(4,1);

            for a = 1:numel(t)
                for b = 1:numel(t)
                    X = t(a);
                    Y = t(b);

                    N = [(1-X)*(1-Y); X*(1-Y); (1-X)*Y; X*Y];

                    dNx = [-(1-Y)/hx; (1-Y)/hx; -Y/hx; Y/hx];
                    dNy = [-(1-X)/hy; -X/hy; (1-X)/hy; X/hy];

                    xq = x0 + hx*X;
                    yq = y0 + hy*Y;

                    Ke = Ke + (dNx*dNx' + dNy*dNy' + c0_fun(xq,yq) * (N*N')) * (hx*hy*wt(a)*wt(b));
                    Fe = Fe + N * f_fun(xq,yq) * (hx*hy*wt(a)*wt(b));
                end
            end

            K(nd,nd) = K(nd,nd) + Ke;
            rhs0(nd) = rhs0(nd) + Fe;
        end
    end

    px = reshape(xg.',[],1);
    py = reshape(yg.',[],1);

    if isempty(BC)
        return;
    end

    u = zeros(Nn,1);

    idx_left = node(ones(my+1,1), (1:(my+1))');
    idx_right = node((mx+1)*ones(my+1,1), (1:(my+1))');
    idx_bottom = node((1:(mx+1))', ones(mx+1,1));
    idx_top = node((1:(mx+1))', (my+1)*ones(mx+1,1));

    idx_c1 = node(1,1);
    idx_c2 = node(mx+1,1);
    idx_c3 = node(1,my+1);
    idx_c4 = node(mx+1,my+1);

    u(idx_left) = BC{1}(px(idx_left), py(idx_left));
    u(idx_right) = BC{2}(px(idx_right), py(idx_right));
    u(idx_bottom) = BC{3}(px(idx_bottom), py(idx_bottom));
    u(idx_top) = BC{4}(px(idx_top), py(idx_top));

    u(idx_c1) = 0.5 * (BC{1}(ax,ay) + BC{3}(ax,ay));
    u(idx_c2) = 0.5 * (BC{2}(bx,ay) + BC{3}(bx,ay));
    u(idx_c3) = 0.5 * (BC{1}(ax,by) + BC{4}(ax,by));
    u(idx_c4) = 0.5 * (BC{2}(bx,by) + BC{4}(bx,by));

    idx_boundary = unique([idx_left; idx_right; idx_bottom; idx_top]);
    idx_interior = (1:Nn).';
    idx_interior(idx_boundary) = [];

    rhs0 = rhs0 - K(:,idx_boundary) * u(idx_boundary);

    u(idx_interior) = K(idx_interior,idx_interior) \ rhs0(idx_interior);
end