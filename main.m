clear;

k=pi;

for div = 4:4
    calculate(2,div,true,@()DtNTest1(k));
end

function calculate(div,divP,plot,test)

[u_ref,f_fun] = test();

N    = 2^divP;

s = cell(N*N,1);
k = 0;
for jy = 1:N
    ay = (jy-1)/N; by = jy/N;
    for ix = 1:N
        ax = (ix-1)/N; bx = ix/N;
        k = k + 1;
        s{k} = Subdomain(div, f_fun, ax, bx, ay, by);
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

[u_global] = reconstructVolumeSolution(s,divP,div,u_ref,nd,u_skel');

[uMono,~,~,pxG,pyG] = get_fem(divP+div,u_ref,f_fun);

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
fprintf("L2 error between discrete solutions: %.5e\n", norm(uMono-u_global));

end
