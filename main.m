clear;

k = 2;

poincareSteklovOperator = "ItI"; 

for div = 3:3
    calculate(4, div, true, @()test1(k), poincareSteklovOperator, k);
end

function calculate(div, divP, plot, test, poincareSteklovOperator, k)

[u_ref,rhs,c0,BC,IBC] = test();

if poincareSteklovOperator == "ItI"
    eta    = k;
    sideBC = IBC;
else
    eta    = 0;
    sideBC = BC;
end

N    = 2^divP;
Xdom = [0, 1];
Ydom = [0, 1];

s = cell(N*N, 1);
cnt = 0;
tol = 1.0e-12;

for jy = 1:N
    ay = Ydom(1) + (jy-1) * (Ydom(2)-Ydom(1)) / N;
    by = Ydom(1) +  jy    * (Ydom(2)-Ydom(1)) / N;
    for ix = 1:N
        ax = Xdom(1) + (ix-1) * (Xdom(2)-Xdom(1)) / N;
        bx = Xdom(1) +  ix    * (Xdom(2)-Xdom(1)) / N;
        cnt = cnt + 1;

        s{cnt} = Subdomain(div, rhs, ax, bx, ay, by, eta, c0, poincareSteklovOperator, @get_sem);

        if abs(s{cnt}.ax - Xdom(1)) < tol
            s{cnt}.setBoundaryCondition(sideBC{1}(s{cnt}.px(s{cnt}.idx(1)), s{cnt}.py(s{cnt}.idx(1))), 1);
        end
        if abs(s{cnt}.bx - Xdom(2)) < tol
            s{cnt}.setBoundaryCondition(sideBC{2}(s{cnt}.px(s{cnt}.idx(2)), s{cnt}.py(s{cnt}.idx(2))), 2);
        end
        if abs(s{cnt}.ay - Ydom(1)) < tol
            s{cnt}.setBoundaryCondition(sideBC{3}(s{cnt}.px(s{cnt}.idx(3)), s{cnt}.py(s{cnt}.idx(3))), 3);
        end
        if abs(s{cnt}.by - Ydom(2)) < tol
            s{cnt}.setBoundaryCondition(sideBC{4}(s{cnt}.px(s{cnt}.idx(4)), s{cnt}.py(s{cnt}.idx(4))), 4);
        end

        ic = s{cnt}.idx(5);
        gc = zeros(numel(ic), 1);
        on_global = (abs(s{cnt}.px(ic) - Xdom(1)) < tol) | (abs(s{cnt}.px(ic) - Xdom(2)) < tol) | ...
                    (abs(s{cnt}.py(ic) - Ydom(1)) < tol) | (abs(s{cnt}.py(ic) - Ydom(2)) < tol);
        if any(on_global)
            if poincareSteklovOperator == "DtN"
                gc(on_global) = BC{1}(s{cnt}.px(ic(on_global)), s{cnt}.py(ic(on_global)));
            else
                cornerExtSides = {[1,3], [2,3], [1,4], [2,4]};
                for c = find(on_global)'
                    for ss = cornerExtSides{c}(:).'
                        sideIsExt = (ss==1 && abs(s{cnt}.ax-Xdom(1))<tol) || ...
                                    (ss==2 && abs(s{cnt}.bx-Xdom(2))<tol) || ...
                                    (ss==3 && abs(s{cnt}.ay-Ydom(1))<tol) || ...
                                    (ss==4 && abs(s{cnt}.by-Ydom(2))<tol);
                        if sideIsExt
                            gc(c) = gc(c) + 0.5 * IBC{ss}(s{cnt}.px(ic(c)), s{cnt}.py(ic(c)));
                        end
                    end
                end
            end
            s{cnt}.setBoundaryCondition(gc, 5);
        end
    end
end

% ================================================================
%  Assemble, solve, reconstruct
% ================================================================


if poincareSteklovOperator == "DtN"
    nd = NestedDissection(divP);
    nd.divide(0);
    nd.calculateReordering(2^div - 1);
    [S,R] = assembleDtN(s, div, nd);
else
    nd = NestedDissectionItI(divP);
    nd.divide(0);
    nd.calculateReordering(2^div - 1);
    [S,R] = assembleItI(s, div, nd, IBC);
end

S2 = S(nd.permutation, nd.permutation);
R2 = R(nd.permutation);

m = 1;
S2inv = Multigrid(S2, nd, length(nd.levels), m);

M   = @(r) S2inv.vcycle(r);
rho = 1 - 1/(2*m + 1)^2;

u_skel = zeros(size(S2, 1), 1);
MR2 = M(R2);
u_skel(nd.permutation) = (1 + 1/rho) * MR2 - (1/rho) * M(S2 * MR2);

fprintf('Residual norm of skeleton solve: %e\n', norm(R2 - S2*u_skel(nd.permutation)));

if poincareSteklovOperator == "DtN"
   % plot_points_in_nd_permutation_order(s, div, nd, S2);
    [u_global, ~, pxG, pyG] = reconstructVolumeSolution(s, divP, div, u_ref, nd, u_skel.');
else
   % plotItIDofIndices(s,div,nd,S2);
    [u_global, ~, pxG, pyG] = reconstructVolumeSolutionItI(s, divP, div, IBC, nd, u_skel);
end

% ================================================================
%  Error reporting
% ================================================================
u_ref_global = u_ref(pxG, pyG);

fprintf("div: %i   divP: %i   (%i x %i subdomains)\n", div, divP, N, N);
fprintf("L2 error u_global vs u_ref (nodal): %.5e\n", norm(u_global - u_ref_global));
fprintf("Rel L2 error u_global vs u_ref:     %.5e\n", norm(u_global - u_ref_global) / norm(u_ref_global));

% ================================================================
%  Compare with monolithic FEM
% ================================================================
[~,~,~,~,~,uMono] = get_fem(divP+div, c0, rhs, Xdom(1), Xdom(2), Ydom(1), Ydom(2), ...
                              2^(divP+div), 2^(divP+div), BC);
fprintf("L2 error between discrete solutions: %.5e\n", norm(uMono - u_global));

% ================================================================
%  Plot
% ================================================================
if plot
    rowsA = N * 2^div + 1;
    X = reshape(pxG, rowsA, rowsA);
    Y = reshape(pyG, rowsA, rowsA);
    Z = reshape(real(u_global), rowsA, rowsA);
    figure;
    surf(X, Y, Z);
    view(3);
end

end

function plotItIDofIndices(s,div,nd,S2)
dofPerFace = 2^div - 1;
nFaces = size(nd.elementsPerFace,1);
totalFaceDofs = 2*nFaces*dofPerFace;

pointDofPerElement = zeros(size(nd.crossPointsPerElement));
pointDofCursor = 0;
for e = 1:size(nd.crossPointsPerElement,1)
    for c = 1:4
        if nd.crossPointsPerElement(e,c) ~= 0
            pointDofCursor = pointDofCursor + 1;
            pointDofPerElement(e,c) = pointDofCursor;
        end
    end
end
totalPointDofs = pointDofCursor;

nTot = totalFaceDofs + totalPointDofs;
xp = nan(nTot,1);
yp = nan(nTot,1);
dx = zeros(nTot,1);
dy = zeros(nTot,1);

h = min(cellfun(@(t) min(t.bx - t.ax, t.by - t.ay), s));
alpha = 0.12;
epsLabel = -0.2*h/(2^div);

for i = 1:nFaces
    for q = 1:2
        b = 2*(i-1) + q;
        ids = (b-1)*dofPerFace + (1:dofPerFace);

        e = nd.elementsPerFace(i,q);
        f = nd.elementSidePerFace(i,q);
        ii = s{e}.idx(f);

        xx = s{e}.px(ii);
        yy = s{e}.py(ii);

        cx = 0.5*(s{e}.ax + s{e}.bx);
        cy = 0.5*(s{e}.ay + s{e}.by);

        xx = (1 - alpha)*xx + alpha*cx;
        yy = (1 - alpha)*yy + alpha*cy;

        xp(ids) = xx;
        yp(ids) = yy;

        dx(ids) = sign(xx - cx)*epsLabel;
        dy(ids) = sign(yy - cy)*epsLabel;
        dx(ids(dx(ids)==0)) = epsLabel;
    end
end

for e = 1:numel(s)
    for c = 1:4
        pd = pointDofPerElement(e,c);
        if pd ~= 0
            id = totalFaceDofs + pd;

            ii = s{e}.idx_corners(c);
            xx = s{e}.px(ii);
            yy = s{e}.py(ii);

            cx = 0.5*(s{e}.ax + s{e}.bx);
            cy = 0.5*(s{e}.ay + s{e}.by);

            xx = (1 - alpha)*xx + alpha*cx;
            yy = (1 - alpha)*yy + alpha*cy;

            xp(id) = xx;
            yp(id) = yy;

            dx(id) = sign(xx - cx)*epsLabel;
            dy(id) = sign(yy - cy)*epsLabel;
            if dx(id) == 0
                dx(id) = epsLabel;
            end
        end
    end
end

p = nd.permutation(:);

figure;
scatter(xp(p),yp(p),30,'filled');
axis equal;
hold on

S2p = S2(1:numel(p),1:numel(p));
[I,J] = find(triu(S2p,1));
for k = 1:numel(I)
    gi = p(I(k));
    gj = p(J(k));
    line([xp(gi) xp(gj)],[yp(gi) yp(gj)], ...
        'Color',[0 0 0], ...
        'LineWidth',0.25);
end

for i = 1:nTot
    text(xp(p(i)) + dx(p(i)), yp(p(i)) + dy(p(i)), sprintf('%d',i), ...
        'FontSize',8, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end

title('Raw ordering');
hold off
end

function plot_points_in_nd_permutation_order_ItI(s, div, nd, varargin)
m = 2^div - 1;
nf = size(nd.elementsPerFace,1);
nic = size(nd.elementsPerCrossPoint,1);
nbc = size(nd.elementsPerBoundaryCrossPoint,1);
nfd = 2*m*nf;
nt = nfd + 4*nic + 2*nbc;

x = nan(nt,1);
y = nan(nt,1);

h = min(cellfun(@(t) min(t.bx-t.ax,t.by-t.ay), s));
e = 0.03*h/(2^div);

for f = 1:nf
    for q = 1:2
        k = (f-1)*2*m + (q-1)*m + (1:m);
        t = s{nd.elementsPerFace(f,q)};
        r = t.idx(nd.elementSidePerFace(f,q));
        x(k) = t.px(r);
        y(k) = t.py(r);
        if nd.elementSidePerFace(f,q) == 1
            x(k) = x(k) + e;
        elseif nd.elementSidePerFace(f,q) == 2
            x(k) = x(k) - e;
        elseif nd.elementSidePerFace(f,q) == 3
            y(k) = y(k) + e;
        else
            y(k) = y(k) - e;
        end
    end
end

for k = 1:nic
    r = nfd + 4*(k-1);
    t = s{nd.elementsPerCrossPoint(k,1)};
    x(r+1) = t.bx + e;
    y(r+1) = t.by + e;
    t = s{nd.elementsPerCrossPoint(k,2)};
    x(r+2) = t.ax - e;
    y(r+2) = t.by + e;
    t = s{nd.elementsPerCrossPoint(k,3)};
    x(r+3) = t.bx + e;
    y(r+3) = t.ay - e;
    t = s{nd.elementsPerCrossPoint(k,4)};
    x(r+4) = t.ax - e;
    y(r+4) = t.ay - e;
end

for k = 1:nbc
    for q = 1:2
        r = nfd + 4*nic + 2*(k-1) + q;
        t = s{nd.elementsPerBoundaryCrossPoint(k,q)};
        c = nd.cornersPerBoundaryCrossPoint(k,q);
        if c == 1
            x(r) = t.ax - e;
            y(r) = t.by + e;
        elseif c == 2
            x(r) = t.bx + e;
            y(r) = t.by + e;
        elseif c == 3
            x(r) = t.bx + e;
            y(r) = t.ay - e;
        else
            x(r) = t.ax - e;
            y(r) = t.ay - e;
        end
    end
end

if nargin > 3
    p = nd.permutation(varargin{1});
else
    p = nd.permutation;
end

figure('Position',[100 100 1400 1400]);
scatter(x(p),y(p),20,'filled');
axis equal;
axis([0 1 0 1]);
hold on

for j = 1:numel(p)
    id = p(j);
    if id <= nfd
        f = floor((id-1)/(2*m)) + 1;
        q = floor(mod(id-1,2*m)/m) + 1;
        side = nd.elementSidePerFace(f,q);
        if side == 1
            dx = 4*e;
            dy = 0;
        elseif side == 2
            dx = -4*e;
            dy = 0;
        elseif side == 3
            dx = 0;
            dy = 4*e;
        else
            dx = 0;
            dy = -4*e;
        end
    elseif id <= nfd + 4*nic
        q = mod(id-nfd-1,4) + 1;
        if q == 1
            dx = 4*e;
            dy = 4*e;
        elseif q == 2
            dx = -4*e;
            dy = 4*e;
        elseif q == 3
            dx = 4*e;
            dy = -4*e;
        else
            dx = -4*e;
            dy = -4*e;
        end
    else
        r = id - nfd - 4*nic;
        k = floor((r-1)/2) + 1;
        q = mod(r-1,2) + 1;
        c = nd.cornersPerBoundaryCrossPoint(k,q);
        if c == 1
            dx = -4*e;
            dy = 4*e;
        elseif c == 2
            dx = 4*e;
            dy = 4*e;
        elseif c == 3
            dx = 4*e;
            dy = -4*e;
        else
            dx = -4*e;
            dy = -4*e;
        end
    end
    text(x(id)+dx,y(id)+dy,sprintf('%d',id),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'FontWeight','bold');
end

hold off
end