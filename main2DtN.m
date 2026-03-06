clear;

k = 2;
tol = 1.0e-14;
[u_exact,rhs,c0,BC,~,~,~,~] = DtNTest1(k);

ay = 0;
by = 1;
axW = 0;
bxW = 1;
axE = 1;
bxE = 2;

for div = 5:5
    sdW = Subdomain(div, rhs, axW, bxW, ay, by, 0, c0, "DtN", @get_sem);
    sdE = Subdomain(div, rhs, axE, bxE, ay, by, 0, c0, "DtN", @get_sem);

    ubW = cell(5,1);
    ubE = cell(5,1);

    ubW{1} = BC{1}(sdW.px(sdW.idx_left), sdW.py(sdW.idx_left));
    ubW{2} = zeros(sdW.b(2),1);
    ubW{3} = BC{3}(sdW.px(sdW.idx_bottom), sdW.py(sdW.idx_bottom));
    ubW{4} = BC{4}(sdW.px(sdW.idx_top), sdW.py(sdW.idx_top));
    ubW{5} = u_exact(sdW.px(sdW.idx_corners), sdW.py(sdW.idx_corners));

    ubE{1} = zeros(sdE.b(1),1);
    ubE{2} = BC{2}(sdE.px(sdE.idx_right), sdE.py(sdE.idx_right));
    ubE{3} = BC{3}(sdE.px(sdE.idx_bottom), sdE.py(sdE.idx_bottom));
    ubE{4} = BC{4}(sdE.px(sdE.idx_top), sdE.py(sdE.idx_top));
    ubE{5} = u_exact(sdE.px(sdE.idx_corners), sdE.py(sdE.idx_corners));

    uI = (sdW.T{2,2} + sdE.T{1,1}) \ ( ...
        (sdW.h{2} - sdW.T{2,1} * ubW{1} - sdW.T{2,3} * ubW{3} - sdW.T{2,4} * ubW{4} - sdW.T{2,5} * ubW{5}) + ...
        (sdE.h{1} - sdE.T{1,2} * ubE{2} - sdE.T{1,3} * ubE{3} - sdE.T{1,4} * ubE{4} - sdE.T{1,5} * ubE{5}) );

    ubW{2} = uI;
    ubE{1} = uI;

    ubWv = [ubW{1}; ubW{2}; ubW{3}; ubW{4}; ubW{5}];
    ubEv = [ubE{1}; ubE{2}; ubE{3}; ubE{4}; ubE{5}];

    uW = zeros(size(sdW.rhs));
    uE = zeros(size(sdE.rhs));

    uW(sdW.idx_interior) = sdW.A \ (sdW.rhs(sdW.idx_interior) - sdW.B * ubWv);
    uE(sdE.idx_interior) = sdE.A \ (sdE.rhs(sdE.idx_interior) - sdE.B * ubEv);

    uW(sdW.idx_boundary) = ubWv;
    uE(sdE.idx_boundary) = ubEv;

    uexW = u_exact(sdW.px, sdW.py);
    uexE = u_exact(sdE.px, sdE.py);

    fprintf("div: %i\n", div);
    fprintf("Rel L2 error = %.15e\n", sqrt(norm(uW-uexW)^2 + norm(uE-uexE)^2) / sqrt(norm(uexW)^2 + norm(uexE)^2));

    nW = 2^div + 1;
    nE = 2^div + 1;

    figure;
    hold on;
    surf(reshape(sdW.px,nW,nW), reshape(sdW.py,nW,nW), reshape(real(uW),nW,nW));
    surf(reshape(sdE.px,nE,nE), reshape(sdE.py,nE,nE), reshape(real(uE),nE,nE));
    view(3);
    hold off;
end