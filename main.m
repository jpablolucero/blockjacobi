clear;

k = 4;

poincareSteklovOperator = "DtN";   % "DtN" or "ItI"

for div = 2:2
    calculate(5, div, true, @()ItITest2(k), poincareSteklovOperator, k);
end

function calculate(div, divP, doPlot, test, poincareSteklovOperator, k)

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

        s{cnt} = Subdomain(div, rhs, ax, bx, ay, by, eta, c0, poincareSteklovOperator, @get_fem);

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

    [S, R] = assembleDtN(s, div, nd);

    S2 = S(nd.permutation, nd.permutation);
    R2 = R(nd.permutation);

    m = 1;
    S2inv = Multigrid(S2, nd, length(nd.levels), m);

    M   = @(r) S2inv.vcycle(r);
    rho = 1 - 1/(2*m + 1)^2;

    u_skel = zeros(size(S2, 1), 1);
    MR2 = M(R2);
    u_skel(nd.permutation) = (1 + 1/rho) * MR2 - (1/rho) * M(S2 * MR2);

    fprintf('Final   L2-norm of GMRES:       %d\n', norm(R2 - S2*u_skel(nd.permutation)));

    [u_global, ~, pxG, pyG] = reconstructVolumeSolution(s, divP, div, u_ref, nd, u_skel.');
else
    [S, R, skel] = assembleItI(s, div, divP, IBC, Xdom, Ydom);

    u_skel = S \ R;

    fprintf('Residual norm of skeleton solve: %e\n', norm(R - S*u_skel));

    [u_global, ~, pxG, pyG] = reconstructVolumeSolutionItI( ...
        s, divP, div, u_ref, u_skel, skel, IBC, eta);
end

% ================================================================
%  Error reporting
% ================================================================
u_ref_global = u_ref(pxG, pyG);

fprintf("div: %i   divP: %i   (%i x %i subdomains)\n", div, divP, N, N);
fprintf("L2 error u_global vs u_ref (nodal): %.5e\n", norm(u_global - u_ref_global));
fprintf("Rel L2 error u_global vs u_ref:     %.5e\n", norm(u_global - u_ref_global) / norm(u_ref_global));

% ================================================================
%  Plot
% ================================================================
if doPlot
    rowsA = N * 2^div + 1;
    X = reshape(pxG, rowsA, rowsA);
    Y = reshape(pyG, rowsA, rowsA);
    Z = reshape(real(u_global), rowsA, rowsA);
    figure;
    surf(X, Y, Z);
    view(3);
end

% ================================================================
%  Compare with monolithic FEM
% ================================================================
[~,~,~,~,~,uMono] = get_fem(divP+div, c0, rhs, Xdom(1), Xdom(2), Ydom(1), Ydom(2), ...
                              2^(divP+div), 2^(divP+div), BC);
fprintf("L2 error between discrete solutions: %.5e\n", norm(uMono - u_global));

end