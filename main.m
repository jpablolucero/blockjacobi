clear;

k=2;

for div = 3:3
    calculate(3,div,true,@()DtNTest1(k));
end

function calculate(div,divP,plot,test);

[u_ref,f_fun,c0,BC,poincareSteklovOperator] = test();

N    = 2^divP;

s = cell(N*N,1);
k = 0;
for jy = 1:N
    ay = (jy-1)/N; by = jy/N;
    for ix = 1:N
        ax = (ix-1)/N; bx = ix/N;
        k = k + 1;
        s{k} = SubdomainSEM(div, f_fun, ax, bx, ay, by, 1, 1, 0, c0, poincareSteklovOperator);
    end 
end

nd = NestedDissection(divP);
nd.divide(0);
nd.calculateReordering(2^div - 1);

[S,R] = assembleDtN(u_ref,s,div,nd);

S2 = S(nd.permutation,nd.permutation);
R2 = R(nd.permutation);

t1 = tic;

m = 1;

S2inv = Multigrid(S2,nd,length(nd.levels),m);

fprintf('Build HPS:\t\t %.6f s\n', toc(t1));
t1 = tic;

M  = @(r) S2inv.vcycle(r);
rho = 1 - 1/(2*m + 1)^2;

u_skel = zeros(size(S2,1),1);

MR2 = M(R2);
u_skel(nd.permutation) = (1 + 1/rho) * MR2 - (1/rho) * M(S2 * MR2);

fprintf('Initial L2-norm of GMRES:       %d\n', norm(M(R2)));
fprintf('Final   L2-norm of GMRES:       %d\n', norm(R2 - S2*u_skel(nd.permutation)));

fprintf('Solve HPS:\t\t %.6f s\n', toc(t1));

% nstart = 1;
% plot_points_in_nd_permutation_order(s, div, nd, S2(nstart:end,nstart:end),nstart:length(R2));

[u_global,u_cells,pxG,pyG] = reconstructVolumeSolution(s,divP,div,u_ref,nd,u_skel');

u_ref_global = u_ref(pxG,pyG);
fprintf("L2 error u_global vs u_ref (nodal): %.5e\n", norm(u_global - u_ref_global));
fprintf("Rel L2 error u_global vs u_ref:     %.5e\n", norm(u_global - u_ref_global)/norm(u_ref_global));

if plot
    rowsA = (2^divP)*(2^div) + 1;
    X = reshape(pxG, rowsA, rowsA);
    Y = reshape(pyG, rowsA, rowsA);
    Z = reshape(u_global, rowsA, rowsA);
    surf(X, Y, Z);
end

[uMono] = get_fem(divP+div,u_ref,f_fun);
fprintf("L2 error between discrete solutions: %.5e\n", norm(uMono-u_global));

end

function plot_points_in_nd_permutation_order(s, div, nd, S2, varargin)
b = 2^div - 1;
m = size(nd.elementsPerFace, 1) * b;
ncp = max(nd.crossPointsPerElement(:));
n = m + ncp;

x = nan(n, 1);
y = nan(n, 1);

for e = 1:numel(s)
    for f = 1:4
        k = nd.facePerElement(e, f);
        if k ~= 0
            ids = (k - 1) * b + (1:b);
            ii = s{e}.idx(f);
            t = isnan(x(ids));
            x(ids(t)) = s{e}.px(ii(t));
            y(ids(t)) = s{e}.py(ii(t));
        end
    end

    for j = 1:4
        kp = nd.crossPointsPerElement(e, j);
        if kp ~= 0
            id = m + kp;
            if isnan(x(id))
                ii = s{e}.idx_corners(j);
                x(id) = s{e}.px(ii);
                y(id) = s{e}.py(ii);
            end
        end
    end
end

if (length(varargin)>0)
    p = nd.permutation(varargin{1});
else
    p = nd.permutation(:);
end

figure;
scatter(x(p), y(p), 20, 'filled');
axis equal;
axis([0 1 0 1]);

% Lines
pp = p(:);                          
S2p = S2(1:numel(pp), 1:numel(pp)); 

[I,J] = find(triu(S2p,1));
hold on
for k = 1:numel(I)
    gi = pp(I(k));  gj = pp(J(k));
    line([x(gi) x(gj)], [y(gi) y(gj)], 'Color',[0 0 0], 'LineWidth',0.25);
end
hold off

for i = 1:numel(p)
    text(x(p(i)), y(p(i)), sprintf('%d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
end