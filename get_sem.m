function [px,py,K,rhs,Mb,u] = get_sem(div,c0_fun,rhs_fun,ax,bx,ay,by,BC)
    if nargin < 8
        BC = [];
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

    if isa(c0_fun,'function_handle')
        c0eval = c0_fun;
    else
        c0const = c0_fun;
        c0eval = @(x,y) c0const + 0.*x;
    end

    p = 2^div;
    quad = Quadrature(p);

    x = quad.Points;
    y = quad.Points;

    xscale = (bx - ax) / 2;
    yscale = (by - ay) / 2;

    xCoords = xscale * x + (bx + ax) / 2;
    yCoords = yscale * y + (by + ay) / 2;

    [X,Y] = meshgrid(xCoords,yCoords);

    px = reshape(X.',[],1);
    py = reshape(Y.',[],1);

    xWeights = xscale * quad.Weights;
    yWeights = yscale * quad.Weights;

    Dx = computeDifferentiationMatrix(x, quad.Polynomials) / xscale;
    Dy = computeDifferentiationMatrix(y, quad.Polynomials) / yscale;

    Mx = spdiags(xWeights,0,numel(xWeights),numel(xWeights));
    My = spdiags(yWeights,0,numel(yWeights),numel(yWeights));

    K = kron(Dx' * Mx * Dx, My) + kron(Mx, Dy' * My * Dy) + c0eval(px,py) .* kron(My,Mx);
    K = full(K);

    rhs = kron(xWeights,yWeights) .* rhs_fun(px,py);

    tol = 1.0e-14;

    idx_left = find(abs(px-ax) < tol & ~(abs(py-ay) < tol | abs(py-by) < tol));
    idx_right = find(abs(px-bx) < tol & ~(abs(py-ay) < tol | abs(py-by) < tol));
    idx_bottom = find(abs(py-ay) < tol & ~(abs(px-ax) < tol | abs(px-bx) < tol));
    idx_top = find(abs(py-by) < tol & ~(abs(px-ax) < tol | abs(px-bx) < tol));
    idx_corners = find(((abs(px-ax) < tol) | (abs(px-bx) < tol)) & ((abs(py-ay) < tol) | (abs(py-by) < tol)));

    Nn = numel(px);

    w = sparse(Nn,1);

    wY = yWeights(2:end-1);
    wX = xWeights(2:end-1);

    w(idx_left) = wY;
    w(idx_right) = wY;
    w(idx_bottom) = wX;
    w(idx_top) = wX;

    ic = idx_corners;

    is_c1 = (abs(px(ic)-ax)<tol) & (abs(py(ic)-ay)<tol);
    is_c2 = (abs(px(ic)-bx)<tol) & (abs(py(ic)-ay)<tol);
    is_c3 = (abs(px(ic)-ax)<tol) & (abs(py(ic)-by)<tol);
    is_c4 = (abs(px(ic)-bx)<tol) & (abs(py(ic)-by)<tol);

    w(ic(is_c1)) = xWeights(1)   + yWeights(1);
    w(ic(is_c2)) = xWeights(end) + yWeights(1);
    w(ic(is_c3)) = xWeights(1)   + yWeights(end);
    w(ic(is_c4)) = xWeights(end) + yWeights(end);

    Mb = spdiags(w,0,Nn,Nn);

    if isempty(BC)
        u = [];
        return;
    end

    u = zeros(Nn,1);

    u(idx_left) = BC{1}(px(idx_left),py(idx_left));
    u(idx_right) = BC{2}(px(idx_right),py(idx_right));
    u(idx_bottom) = BC{3}(px(idx_bottom),py(idx_bottom));
    u(idx_top) = BC{4}(px(idx_top),py(idx_top));

    for j = 1:numel(idx_corners)
        id = idx_corners(j);

        vals = [];

        if abs(px(id)-ax) < tol
            vals(end+1,1) = BC{1}(px(id),py(id));
        end
        if abs(px(id)-bx) < tol
            vals(end+1,1) = BC{2}(px(id),py(id));
        end
        if abs(py(id)-ay) < tol
            vals(end+1,1) = BC{3}(px(id),py(id));
        end
        if abs(py(id)-by) < tol
            vals(end+1,1) = BC{4}(px(id),py(id));
        end

        u(id) = mean(vals);
    end

    idx_boundary = [idx_left; idx_right; idx_bottom; idx_top; idx_corners];
    idx_interior = (1:Nn).';
    idx_interior(idx_boundary) = [];

    rhs = rhs - K(:,idx_boundary) * u(idx_boundary);

    u(idx_interior) = K(idx_interior,idx_interior) \ rhs(idx_interior);
end

function D = computeDifferentiationMatrix(x, P)
    N1 = length(x);

    X = repmat(x, 1, N1);
    Xdiff = X - X' + eye(N1);

    L = repmat(P(:, end), 1, N1);
    L(1:(N1+1):N1*N1) = 1;

    D = (L ./ (Xdiff .* L'));
    D(1:(N1+1):N1*N1) = 0;

    N = N1 - 1;
    D(1, 1) = -(N1 * N) / 4;
    D(end, end) = (N1 * N) / 4;
end