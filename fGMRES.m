function [x,iter,resids] = fGMRES(A, b, tol, varargin)
% Flexible GMRES method
%
%   [x,iter,resids] = fGMRES(A, b, tol, varargin)
%
% A may be
%   - a matrix
%   - a function handle @(x,tolA) A*x
%
% P may be
%   - empty
%   - a matrix
%   - a function handle @(x,tolP) P*x
%
% Parameters:
%   'relaxation' [0]
%   'max_iters'  [2]   number of restart cycles
%   'restart'    [20]  inner iterations per cycle
%   'x0'         [zeros]
%   'verb'       [2]
%   'tol_exit'   [tol]
%   'P'          []
%
% Returns:
%   x      approximate solution
%   iter   [outer_cycle, inner_iter]
%   resids residual history per inner step, size restart x outer_cycles

relaxation = 0;
max_iters = 2;
restart = 20;
x = [];
verb = 2;
tol_exit = tol;
P = [];

for i = 1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'relaxation'
            relaxation = varargin{i+1};
        case 'max_iters'
            max_iters = varargin{i+1};
        case 'restart'
            restart = varargin{i+1};
        case 'x0'
            x = varargin{i+1};
        case 'verb'
            verb = varargin{i+1};
        case 'tol_exit'
            tol_exit = varargin{i+1};
        case 'p'
            P = varargin{i+1};
        otherwise
            error('Unknown tuning parameter "%s"', varargin{i});
    end
end

if isa(A, 'double')
    A0 = A;
    A = @(x,tolA) A0 * x;
end

if ~isempty(P)
    if isa(P, 'double')
        P0 = P;
        P = @(x,tolP) P0 * x;
    end
end

n = size(b,1);

if restart > n
    restart = n;
end

beta0 = norm(b);

if beta0 == 0
    x = zeros(n,1);
    iter = [0, 0];
    resids = 0;
    return
end

if isempty(x)
    x = zeros(n,1);
end

if nargout > 2
    resids = zeros(restart, max_iters);
end

t_gmres_start = tic;

outer_done = 0;
inner_done = 0;
last_resid = norm(b - A(x, tol)) / beta0;

for it = 1:max_iters
    outer_done = it;

    r = b - A(x, tol);
    beta = norm(r);
    resid0 = beta / beta0;

    if verb > 1
        fprintf('==fgmres== iter=%d, |r|=%g, relres=%3.3e, time=%g\n', ...
            it-1, beta, resid0, toc(t_gmres_start));
    end

    if resid0 < tol_exit
        inner_done = 0;
        last_resid = resid0;
        break
    end

    V = cell(restart + 1, 1);
    Z = cell(restart, 1);
    H = zeros(restart + 1, restart);
    cs = zeros(restart, 1);
    sn = zeros(restart, 1);
    g = zeros(restart + 1, 1);

    V{1} = r / beta;
    g(1) = beta;

    happy_breakdown = false;

    for j = 1:restart
        inner_done = j;

        if beta == 0
            tol_kryl = tol;
        else
            err = abs(g(j)) / beta;
            err = max(err, eps);
            tol_kryl = tol / (err^relaxation);
        end

        if isempty(P)
            Z{j} = V{j};
        else
            Z{j} = P(V{j}, tol_kryl);
        end

        w = A(Z{j}, tol_kryl);

        for i = 1:j
            H(i,j) = V{i}' * w;
            w = w - H(i,j) * V{i};
        end

        H(j+1,j) = norm(w);

        if H(j+1,j) <= eps * norm(H(1:j,j))
            happy_breakdown = true;
            H(j+1,j) = 0;
        else
            V{j+1} = w / H(j+1,j);
        end

        for i = 1:j-1
            tmp = cs(i) * H(i,j) + sn(i) * H(i+1,j);
            H(i+1,j) = -conj(sn(i)) * H(i,j) + cs(i) * H(i+1,j);
            H(i,j) = tmp;
        end

        [cs(j), sn(j)] = complex_givens(H(j,j), H(j+1,j));

        H(j,j) = cs(j) * H(j,j) + sn(j) * H(j+1,j);
        H(j+1,j) = 0;

        tmp = cs(j) * g(j) + sn(j) * g(j+1);
        g(j+1) = -conj(sn(j)) * g(j) + cs(j) * g(j+1);
        g(j) = tmp;

        resid = abs(g(j+1)) / beta0;
        last_resid = resid;

        if nargout > 2
            resids(j,it) = resid;
        end

        if verb > 1
            fprintf('==fgmres== iter=[%d,%d], resid=%3.3e, time=%g\n', ...
                it, j, resid, toc(t_gmres_start));
        end

        if resid < tol_exit
            break
        end

        if happy_breakdown
            break
        end
    end

    y = H(1:inner_done,1:inner_done) \ g(1:inner_done);

    dx = zeros(n,1);
    for i = 1:inner_done
        dx = dx + Z{i} * y(i);
    end

    x = x + dx;

    r = b - A(x, tol);
    last_resid = norm(r) / beta0;

    if nargout > 2
        resids(inner_done:restart,it) = last_resid;
    end

    if verb > 0
        fprintf('==fgmres== cycle=%d done, inner=%d, resid=%3.3e, time=%g\n', ...
            it, inner_done, last_resid, toc(t_gmres_start));
    end

    if last_resid < tol_exit
        break
    end
end

iter = [outer_done, inner_done];

if nargout > 2
    resids = resids(:,1:outer_done);
end

if verb > 0
    fprintf('==fgmres== iters: %d outer + %d inner, resid: %3.3e, time: %g\n', ...
        max(outer_done-1,0), inner_done, last_resid, toc(t_gmres_start));
end

end

function [c,s] = complex_givens(a, b)
if b == 0
    c = 1;
    s = 0;
    return
end

if a == 0
    c = 0;
    s = 1;
    return
end

scale = abs(a) + abs(b);
normab = scale * sqrt(abs(a / scale)^2 + abs(b / scale)^2);
alpha = a / abs(a);

c = abs(a) / normab;
s = alpha * conj(b) / normab;
end