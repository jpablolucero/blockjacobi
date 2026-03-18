clear;

k = 2;
eta = k;
[u_exact,rhs,c0,BC,IBC] = DtNTest1(k);
% [u_exact,rhs,c0,BC,IBC] = ItITest1(k);

ay = 0; by = 1;
axW = 0; bxW = 1;
axE = 1; bxE = 2;

for div = 5:5
    sW = Subdomain(div, rhs, axW, bxW, ay, by, 0, c0, "DtN", @get_fem);
    sE = Subdomain(div, rhs, axE, bxE, ay, by, 0, c0, "DtN", @get_fem);

    allEdgesSizeW = sum(sW.b(1:4));
    allEdgesSizeE = sum(sE.b(1:4));

    DtNHW = cell2mat(sW.T);
    DtNHE = cell2mat(sE.T);
    MbW   = full(sW.Mb);
    MbE   = full(sE.Mb);

    ItIHW = MbW \ (DtNHW - 1i*eta*MbW) * ((DtNHW + 1i*eta*MbW) \ MbW);
    ItIPW = -(2i*eta) * ((DtNHW + 1i*eta*MbW) \ cell2mat(sW.h));
    ItIHE = MbE \ (DtNHE - 1i*eta*MbE) * ((DtNHE + 1i*eta*MbE) \ MbE);
    ItIPE = -(2i*eta) * ((DtNHE + 1i*eta*MbE) \ cell2mat(sE.h));

    pxWc = sW.px(sW.idx_corners);  pyWc = sW.py(sW.idx_corners);
    pxEc = sE.px(sE.idx_corners);  pyEc = sE.py(sE.idx_corners);

    iiW = {IBC{1}(sW.px(sW.idx_left),   sW.py(sW.idx_left));
           zeros(sW.b(2), 1);
           IBC{3}(sW.px(sW.idx_bottom), sW.py(sW.idx_bottom));
           IBC{4}(sW.px(sW.idx_top),    sW.py(sW.idx_top));
           [0.5*(IBC{1}(pxWc(1),pyWc(1)) + IBC{3}(pxWc(1),pyWc(1)));
            0.5*IBC{3}(pxWc(2),pyWc(2));
            0.5*(IBC{1}(pxWc(3),pyWc(3)) + IBC{4}(pxWc(3),pyWc(3)));
            0.5*IBC{4}(pxWc(4),pyWc(4))]};

    iiE = {zeros(sE.b(1), 1);
           IBC{2}(sE.px(sE.idx_right),  sE.py(sE.idx_right));
           IBC{3}(sE.px(sE.idx_bottom), sE.py(sE.idx_bottom));
           IBC{4}(sE.px(sE.idx_top),    sE.py(sE.idx_top));
           [0.5*IBC{3}(pxEc(1),pyEc(1));
            0.5*(IBC{2}(pxEc(2),pyEc(2)) + IBC{3}(pxEc(2),pyEc(2)));
            0.5*IBC{4}(pxEc(3),pyEc(3));
            0.5*(IBC{2}(pxEc(4),pyEc(4)) + IBC{4}(pxEc(4),pyEc(4)))]};

    m   = sW.b(2);
    N_i = m + 2;
    N_W = sum(sW.b);
    N_E = sum(sE.b);

    iIface_W = [sW.b(1)+(1:m)'; allEdgesSizeW+2; allEdgesSizeW+4];
    iIface_E = [(1:sE.b(1))';   allEdgesSizeE+1; allEdgesSizeE+3];

    B_W = sparse(N_W, N_i);
    B_W(iIface_W(1:m), 1:m) = speye(m);
    B_W(allEdgesSizeW+2, m+1) = 0.5;
    B_W(allEdgesSizeW+4, m+2) = 0.5;

    B_E = sparse(N_E, N_i);
    B_E(iIface_E(1:m), 1:m) = speye(m);
    B_E(allEdgesSizeE+1, m+1) = 0.5;
    B_E(allEdgesSizeE+3, m+2) = 0.5;

    F_W = sparse(N_i, N_i);  F_W(m+1,m+1) = 0.5;  F_W(m+2,m+2) = 0.5;
    F_E = sparse(N_i, N_i);  F_E(m+1,m+1) = 0.5;  F_E(m+2,m+2) = 0.5;

    g_W = zeros(N_i,1);
    g_W(m+1) = -0.5*IBC{3}(pxWc(2),pyWc(2));
    g_W(m+2) = -0.5*IBC{4}(pxWc(4),pyWc(4));

    g_E = zeros(N_i,1);
    g_E(m+1) = -0.5*IBC{3}(pxEc(1),pyEc(1));
    g_E(m+2) = -0.5*IBC{4}(pxEc(3),pyEc(3));

    iiWv = cell2mat(iiW);
    iiEv = cell2mat(iiE);

    T_eff_W = ItIHW(iIface_W,:)*B_W + F_W;
    h_eff_W = ItIHW(iIface_W,:)*iiWv + ItIPW(iIface_W) + g_W;

    T_eff_E = ItIHE(iIface_E,:)*B_E + F_E;
    h_eff_E = ItIHE(iIface_E,:)*iiEv + ItIPE(iIface_E) + g_E;

    oI = [speye(N_i), T_eff_W; T_eff_E, speye(N_i)] \ [h_eff_W; h_eff_E];

    iiWv = iiWv + B_W*(-oI(N_i+1:end));
    iiEv = iiEv + B_E*(-oI(1:N_i));

    ooWv = ItIHW*iiWv + ItIPW;
    ooEv = ItIHE*iiEv + ItIPE;

    ubWv = (iiWv - ooWv)/(2i*eta);
    ubEv = (iiEv - ooEv)/(2i*eta);

    uW = zeros(size(sW.rhs));
    uE = zeros(size(sE.rhs));

    uW(sW.idx_interior) = sW.A \ (sW.rhs(sW.idx_interior) - sW.B * ubWv);
    uE(sE.idx_interior) = sE.A \ (sE.rhs(sE.idx_interior) - sE.B * ubEv);

    uW(sW.idx_boundary) = ubWv;
    uE(sE.idx_boundary) = ubEv;

    uexW = u_exact(sW.px, sW.py);
    uexE = u_exact(sE.px, sE.py);

    fprintf("div: %i\n", div);
    fprintf("Rel L2 error = %.15e\n", sqrt(norm(uW-uexW)^2 + norm(uE-uexE)^2) / sqrt(norm(uexW)^2 + norm(uexE)^2));

    nW = 2^div + 1;
    nE = 2^div + 1;

    figure; hold on;
    surf(reshape(sW.px,nW,nW), reshape(sW.py,nW,nW), reshape(real(uW),nW,nW));
    surf(reshape(sE.px,nE,nE), reshape(sE.py,nE,nE), reshape(real(uE),nE,nE));
    view(3); hold off;
end