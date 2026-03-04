clear;

k = 2;
eta = k;
tol = 1.0e-14;
[u_exact,rhs,c0,BC,~,IBC,~,~] = DtNTest1(k);

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

    iiW = cell(5,1);
    iiE = cell(5,1);

    iiW{1} = IBC{1}(sdW.px(sdW.idx_left),    sdW.py(sdW.idx_left));
    iiW{2} = zeros(sdW.b(2),1);
    iiW{3} = IBC{3}(sdW.px(sdW.idx_bottom),  sdW.py(sdW.idx_bottom));
    iiW{4} = IBC{4}(sdW.px(sdW.idx_top),     sdW.py(sdW.idx_top));

    iiE{1} = zeros(sdE.b(1),1);
    iiE{2} = IBC{2}(sdE.px(sdE.idx_right),   sdE.py(sdE.idx_right));
    iiE{3} = IBC{3}(sdE.px(sdE.idx_bottom),  sdE.py(sdE.idx_bottom));
    iiE{4} = IBC{4}(sdE.px(sdE.idx_top),     sdE.py(sdE.idx_top));

    TWfull = sdW.Mb\cell2mat(sdW.T(1:5,1:5));
    TEfull = sdW.Mb\cell2mat(sdE.T(1:5,1:5));

    hWfull = sdW.Mb\cell2mat(sdW.h(1:5));
    hEfull = sdW.Mb\cell2mat(sdE.h(1:5));

    nbW4 = sum(sdW.b(1:4));
    nbE4 = sum(sdE.b(1:4));

    gW5 = BC{5}(sdW.px(sdW.idx_corners), sdW.py(sdW.idx_corners));
    gE5 = BC{5}(sdE.px(sdE.idx_corners), sdE.py(sdE.idx_corners));

    gW5([1,3]) = 0;
    gE5([2,4]) = 0;

    gW = [zeros(nbW4,1); gW5];
    gE = [zeros(nbE4,1); gE5];

    hWfull = hWfull - TWfull*gW;
    hEfull = hEfull - TEfull*gE;

    keepW = [(1:nbW4), nbW4+1, nbW4+3].';
    keepE = [(1:nbE4), nbE4+2, nbE4+4].';

    TW = TWfull(keepW,keepW);
    TE = TEfull(keepE,keepE);

    hW = hWfull(keepW);
    hE = hEfull(keepE);

    bW = [sdW.b(1:4), 2];
    bE = [sdE.b(1:4), 2];

    SWv = (TW - 1i*eta*speye(size(TW,1))) / (TW + 1i*eta*speye(size(TW,1)));
    SEv = (TE - 1i*eta*speye(size(TE,1))) / (TE + 1i*eta*speye(size(TE,1)));

    cWv = -(2i*eta) * ((TW + 1i*eta*speye(size(TW,1))) \ hW);
    cEv = -(2i*eta) * ((TE + 1i*eta*speye(size(TE,1))) \ hE);

    SW = mat2cell(SWv, bW.', bW.');  cW = mat2cell(cWv, bW.', 1);
    SE = mat2cell(SEv, bE.', bE.');  cE = mat2cell(cEv, bE.', 1);

    xW1 = sdW.px(sdW.idx_corners(1));  yW1 = sdW.py(sdW.idx_corners(1));
    xW3 = sdW.px(sdW.idx_corners(3));  yW3 = sdW.py(sdW.idx_corners(3));
    iiW{5} = [0.5*(IBC{1}(xW1,yW1) + IBC{3}(xW1,yW1));
              0.5*(IBC{1}(xW3,yW3) + IBC{4}(xW3,yW3))];

    xE2 = sdE.px(sdE.idx_corners(2));  yE2 = sdE.py(sdE.idx_corners(2));
    xE4 = sdE.px(sdE.idx_corners(4));  yE4 = sdE.py(sdE.idx_corners(4));
    iiE{5} = [0.5*(IBC{2}(xE2,yE2) + IBC{3}(xE2,yE2));
              0.5*(IBC{2}(xE4,yE4) + IBC{4}(xE4,yE4))];

    rW = SW{2,1}*iiW{1} + SW{2,3}*iiW{3} + SW{2,4}*iiW{4} + SW{2,5}*iiW{5} + cW{2};
    rE = SE{1,2}*iiE{2} + SE{1,3}*iiE{3} + SE{1,4}*iiE{4} + SE{1,5}*iiE{5} + cE{1};

    m  = sdW.b(2);
    oI = [speye(m), SW{2,2};
          SE{1,1},  speye(m)] \ [rW; rE];

    iiW{2} = -oI(m+(1:m));
    iiE{1} = -oI(1:m);

    ooW = cell(5,1);
    ooE = cell(5,1);

    ooW{1} = SW{1,1}*iiW{1} + SW{1,2}*iiW{2} + SW{1,3}*iiW{3} + SW{1,4}*iiW{4} + SW{1,5}*iiW{5} + cW{1};
    ooW{2} = SW{2,1}*iiW{1} + SW{2,2}*iiW{2} + SW{2,3}*iiW{3} + SW{2,4}*iiW{4} + SW{2,5}*iiW{5} + cW{2};
    ooW{3} = SW{3,1}*iiW{1} + SW{3,2}*iiW{2} + SW{3,3}*iiW{3} + SW{3,4}*iiW{4} + SW{3,5}*iiW{5} + cW{3};
    ooW{4} = SW{4,1}*iiW{1} + SW{4,2}*iiW{2} + SW{4,3}*iiW{3} + SW{4,4}*iiW{4} + SW{4,5}*iiW{5} + cW{4};
    ooW{5} = SW{5,1}*iiW{1} + SW{5,2}*iiW{2} + SW{5,3}*iiW{3} + SW{5,4}*iiW{4} + SW{5,5}*iiW{5} + cW{5};

    ooE{1} = SE{1,1}*iiE{1} + SE{1,2}*iiE{2} + SE{1,3}*iiE{3} + SE{1,4}*iiE{4} + SE{1,5}*iiE{5} + cE{1};
    ooE{2} = SE{2,1}*iiE{1} + SE{2,2}*iiE{2} + SE{2,3}*iiE{3} + SE{2,4}*iiE{4} + SE{2,5}*iiE{5} + cE{2};
    ooE{3} = SE{3,1}*iiE{1} + SE{3,2}*iiE{2} + SE{3,3}*iiE{3} + SE{3,4}*iiE{4} + SE{3,5}*iiE{5} + cE{3};
    ooE{4} = SE{4,1}*iiE{1} + SE{4,2}*iiE{2} + SE{4,3}*iiE{3} + SE{4,4}*iiE{4} + SE{4,5}*iiE{5} + cE{4};
    ooE{5} = SE{5,1}*iiE{1} + SE{5,2}*iiE{2} + SE{5,3}*iiE{3} + SE{5,4}*iiE{4} + SE{5,5}*iiE{5} + cE{5};

    ubW{1} = (iiW{1} - ooW{1})/(2i*eta);
    ubW{2} = (iiW{2} - ooW{2})/(2i*eta);
    ubW{3} = (iiW{3} - ooW{3})/(2i*eta);
    ubW{4} = (iiW{4} - ooW{4})/(2i*eta);

    ubE{1} = (iiE{1} - ooE{1})/(2i*eta);
    ubE{2} = (iiE{2} - ooE{2})/(2i*eta);
    ubE{3} = (iiE{3} - ooE{3})/(2i*eta);
    ubE{4} = (iiE{4} - ooE{4})/(2i*eta);

    ubW5 = BC{5}(sdW.px(sdW.idx_corners), sdW.py(sdW.idx_corners));
    ubE5 = BC{5}(sdE.px(sdE.idx_corners), sdE.py(sdE.idx_corners));

    ubW5([1,3]) = (iiW{5} - ooW{5})/(2i*eta);
    ubE5([2,4]) = (iiE{5} - ooE{5})/(2i*eta);

    ubW{5} = ubW5;
    ubE{5} = ubE5;

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