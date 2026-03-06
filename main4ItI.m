clear;

k = 2;
eta = k;
[u_exact, rhs, c0, ~, ~, IBC, ~, ~] = DtNTest1(k);

for div = 5:5
    sSW = Subdomain(div, rhs, 0.0, 0.5, 0.0, 0.5, 0, c0, "DtN", @get_sem);
    sSE = Subdomain(div, rhs, 0.5, 1.0, 0.0, 0.5, 0, c0, "DtN", @get_sem);
    sNW = Subdomain(div, rhs, 0.0, 0.5, 0.5, 1.0, 0, c0, "DtN", @get_sem);
    sNE = Subdomain(div, rhs, 0.5, 1.0, 0.5, 1.0, 0, c0, "DtN", @get_sem);

    HdtNSW = cell2mat(sSW.T);
    HdtNSE = cell2mat(sSE.T);
    HdtNNW = cell2mat(sNW.T);
    HdtNNE = cell2mat(sNE.T);

    MbSW = full(sSW.Mb);
    MbSE = full(sSE.Mb);
    MbNW = full(sNW.Mb);
    MbNE = full(sNE.Mb);

    PSW = -(2i * eta) * ((HdtNSW + 1i * eta * MbSW) \ cell2mat(sSW.h));
    PSE = -(2i * eta) * ((HdtNSE + 1i * eta * MbSE) \ cell2mat(sSE.h));
    PNW = -(2i * eta) * ((HdtNNW + 1i * eta * MbNW) \ cell2mat(sNW.h));
    PNE = -(2i * eta) * ((HdtNNE + 1i * eta * MbNE) \ cell2mat(sNE.h));

    HSW = MbSW \ (HdtNSW - 1i * eta * MbSW) * ((HdtNSW + 1i * eta * MbSW) \ MbSW);
    HSE = MbSE \ (HdtNSE - 1i * eta * MbSE) * ((HdtNSE + 1i * eta * MbSE) \ MbSE);
    HNW = MbNW \ (HdtNNW - 1i * eta * MbNW) * ((HdtNNW + 1i * eta * MbNW) \ MbNW);
    HNE = MbNE \ (HdtNNE - 1i * eta * MbNE) * ((HdtNNE + 1i * eta * MbNE) \ MbNE);

    tol = 1.0e-14;

    mSW = sSW.b(1);
    nSW = 2 * mSW + 3;
    NSW = sum(sSW.b);
    ESW = sum(sSW.b(1:4));
    icSW = sSW.idx_corners;
    pxcSW = sSW.px(icSW);
    pycSW = sSW.py(icSW);

    jc1SW = find((abs(pxcSW - sSW.ax) < tol) & (abs(pycSW - sSW.ay) < tol));
    jc2SW = find((abs(pxcSW - sSW.bx) < tol) & (abs(pycSW - sSW.ay) < tol));
    jc3SW = find((abs(pxcSW - sSW.ax) < tol) & (abs(pycSW - sSW.by) < tol));
    jc4SW = find((abs(pxcSW - sSW.bx) < tol) & (abs(pycSW - sSW.by) < tol));

    r1SW = ESW + jc1SW;
    r2SW = ESW + jc2SW;
    r3SW = ESW + jc3SW;
    r4SW = ESW + jc4SW;

    ii0SW = zeros(NSW, 1);
    BSW = sparse(NSW, nSW);
    FSW = sparse(nSW, nSW);
    gSW = zeros(nSW, 1);

    ii0SW(1:mSW) = IBC{1}(sSW.px(sSW.idx_left), sSW.py(sSW.idx_left));
    ii0SW(2 * mSW + (1:mSW)) = IBC{3}(sSW.px(sSW.idx_bottom), sSW.py(sSW.idx_bottom));
    ii0SW(r1SW) = 0.5 * (IBC{1}(sSW.px(icSW(jc1SW)), sSW.py(icSW(jc1SW))) + IBC{3}(sSW.px(icSW(jc1SW)), sSW.py(icSW(jc1SW))));
    ii0SW(r2SW) = 0.5 * IBC{3}(sSW.px(icSW(jc2SW)), sSW.py(icSW(jc2SW)));
    ii0SW(r3SW) = 0.5 * IBC{1}(sSW.px(icSW(jc3SW)), sSW.py(icSW(jc3SW)));

    BSW(mSW + (1:mSW), 1:mSW) = speye(mSW);
    BSW(3 * mSW + (1:mSW), mSW + (1:mSW)) = speye(mSW);
    BSW(r2SW, 2 * mSW + 1) = 0.5;
    BSW(r3SW, 2 * mSW + 2) = 0.5;
    BSW(r4SW, 2 * mSW + 3) = 1.0;

    FSW(2 * mSW + 1, 2 * mSW + 1) = 0.5;
    FSW(2 * mSW + 2, 2 * mSW + 2) = 0.5;
    gSW(2 * mSW + 1) = -0.5 * IBC{3}(sSW.px(icSW(jc2SW)), sSW.py(icSW(jc2SW)));
    gSW(2 * mSW + 2) = -0.5 * IBC{1}(sSW.px(icSW(jc3SW)), sSW.py(icSW(jc3SW)));

    rowsSW = [mSW + (1:mSW), 3 * mSW + (1:mSW), r2SW, r3SW, r4SW];

    LSW.T = HSW(rowsSW, :) * BSW + FSW;
    LSW.h = HSW(rowsSW, :) * ii0SW + PSW(rowsSW) + gSW;
    LSW.B = BSW;
    LSW.ii0 = ii0SW;
    LSW.n = nSW;
    LSW.m = mSW;
    LSW.mc = full(sSW.Mb(rowsSW(end), rowsSW(end)));
    LSW.iR = 1:mSW;
    LSW.iT = mSW + (1:mSW);
    LSW.iBR = 2 * mSW + 1;
    LSW.iTL = 2 * mSW + 2;
    LSW.iC = 2 * mSW + 3;

    mSE = sSE.b(1);
    nSE = 2 * mSE + 3;
    NSE = sum(sSE.b);
    ESE = sum(sSE.b(1:4));
    icSE = sSE.idx_corners;
    pxcSE = sSE.px(icSE);
    pycSE = sSE.py(icSE);

    jc1SE = find((abs(pxcSE - sSE.ax) < tol) & (abs(pycSE - sSE.ay) < tol));
    jc2SE = find((abs(pxcSE - sSE.bx) < tol) & (abs(pycSE - sSE.ay) < tol));
    jc3SE = find((abs(pxcSE - sSE.ax) < tol) & (abs(pycSE - sSE.by) < tol));
    jc4SE = find((abs(pxcSE - sSE.bx) < tol) & (abs(pycSE - sSE.by) < tol));

    r1SE = ESE + jc1SE;
    r2SE = ESE + jc2SE;
    r3SE = ESE + jc3SE;
    r4SE = ESE + jc4SE;

    ii0SE = zeros(NSE, 1);
    BSE = sparse(NSE, nSE);
    FSE = sparse(nSE, nSE);
    gSE = zeros(nSE, 1);

    ii0SE(mSE + (1:mSE)) = IBC{2}(sSE.px(sSE.idx_right), sSE.py(sSE.idx_right));
    ii0SE(2 * mSE + (1:mSE)) = IBC{3}(sSE.px(sSE.idx_bottom), sSE.py(sSE.idx_bottom));
    ii0SE(r1SE) = 0.5 * IBC{3}(sSE.px(icSE(jc1SE)), sSE.py(icSE(jc1SE)));
    ii0SE(r2SE) = 0.5 * (IBC{2}(sSE.px(icSE(jc2SE)), sSE.py(icSE(jc2SE))) + IBC{3}(sSE.px(icSE(jc2SE)), sSE.py(icSE(jc2SE))));
    ii0SE(r4SE) = 0.5 * IBC{2}(sSE.px(icSE(jc4SE)), sSE.py(icSE(jc4SE)));

    BSE(1:mSE, 1:mSE) = speye(mSE);
    BSE(3 * mSE + (1:mSE), mSE + (1:mSE)) = speye(mSE);
    BSE(r1SE, 2 * mSE + 1) = 0.5;
    BSE(r4SE, 2 * mSE + 2) = 0.5;
    BSE(r3SE, 2 * mSE + 3) = 1.0;

    FSE(2 * mSE + 1, 2 * mSE + 1) = 0.5;
    FSE(2 * mSE + 2, 2 * mSE + 2) = 0.5;
    gSE(2 * mSE + 1) = -0.5 * IBC{3}(sSE.px(icSE(jc1SE)), sSE.py(icSE(jc1SE)));
    gSE(2 * mSE + 2) = -0.5 * IBC{2}(sSE.px(icSE(jc4SE)), sSE.py(icSE(jc4SE)));

    rowsSE = [1:mSE, 3 * mSE + (1:mSE), r1SE, r4SE, r3SE];

    LSE.T = HSE(rowsSE, :) * BSE + FSE;
    LSE.h = HSE(rowsSE, :) * ii0SE + PSE(rowsSE) + gSE;
    LSE.B = BSE;
    LSE.ii0 = ii0SE;
    LSE.n = nSE;
    LSE.m = mSE;
    LSE.mc = full(sSE.Mb(rowsSE(end), rowsSE(end)));
    LSE.iT = mSE + (1:mSE);
    LSE.iBM = 2 * mSE + 1;
    LSE.iRM = 2 * mSE + 2;
    LSE.iC = 2 * mSE + 3;

    mNW = sNW.b(1);
    nNW = 2 * mNW + 3;
    NNW = sum(sNW.b);
    ENW = sum(sNW.b(1:4));
    icNW = sNW.idx_corners;
    pxcNW = sNW.px(icNW);
    pycNW = sNW.py(icNW);

    jc1NW = find((abs(pxcNW - sNW.ax) < tol) & (abs(pycNW - sNW.ay) < tol));
    jc2NW = find((abs(pxcNW - sNW.bx) < tol) & (abs(pycNW - sNW.ay) < tol));
    jc3NW = find((abs(pxcNW - sNW.ax) < tol) & (abs(pycNW - sNW.by) < tol));
    jc4NW = find((abs(pxcNW - sNW.bx) < tol) & (abs(pycNW - sNW.by) < tol));

    r1NW = ENW + jc1NW;
    r2NW = ENW + jc2NW;
    r3NW = ENW + jc3NW;
    r4NW = ENW + jc4NW;

    ii0NW = zeros(NNW, 1);
    BNW = sparse(NNW, nNW);
    FNW = sparse(nNW, nNW);
    gNW = zeros(nNW, 1);

    ii0NW(1:mNW) = IBC{1}(sNW.px(sNW.idx_left), sNW.py(sNW.idx_left));
    ii0NW(3 * mNW + (1:mNW)) = IBC{4}(sNW.px(sNW.idx_top), sNW.py(sNW.idx_top));
    ii0NW(r1NW) = 0.5 * IBC{1}(sNW.px(icNW(jc1NW)), sNW.py(icNW(jc1NW)));
    ii0NW(r3NW) = 0.5 * (IBC{1}(sNW.px(icNW(jc3NW)), sNW.py(icNW(jc3NW))) + IBC{4}(sNW.px(icNW(jc3NW)), sNW.py(icNW(jc3NW))));
    ii0NW(r4NW) = 0.5 * IBC{4}(sNW.px(icNW(jc4NW)), sNW.py(icNW(jc4NW)));

    BNW(mNW + (1:mNW), 1:mNW) = speye(mNW);
    BNW(2 * mNW + (1:mNW), mNW + (1:mNW)) = speye(mNW);
    BNW(r1NW, 2 * mNW + 1) = 0.5;
    BNW(r4NW, 2 * mNW + 2) = 0.5;
    BNW(r2NW, 2 * mNW + 3) = 1.0;

    FNW(2 * mNW + 1, 2 * mNW + 1) = 0.5;
    FNW(2 * mNW + 2, 2 * mNW + 2) = 0.5;
    gNW(2 * mNW + 1) = -0.5 * IBC{1}(sNW.px(icNW(jc1NW)), sNW.py(icNW(jc1NW)));
    gNW(2 * mNW + 2) = -0.5 * IBC{4}(sNW.px(icNW(jc4NW)), sNW.py(icNW(jc4NW)));

    rowsNW = [mNW + (1:mNW), 2 * mNW + (1:mNW), r1NW, r4NW, r2NW];

    LNW.T = HNW(rowsNW, :) * BNW + FNW;
    LNW.h = HNW(rowsNW, :) * ii0NW + PNW(rowsNW) + gNW;
    LNW.B = BNW;
    LNW.ii0 = ii0NW;
    LNW.n = nNW;
    LNW.m = mNW;
    LNW.mc = full(sNW.Mb(rowsNW(end), rowsNW(end)));
    LNW.iR = 1:mNW;
    LNW.iB = mNW + (1:mNW);
    LNW.iLM = 2 * mNW + 1;
    LNW.iTM = 2 * mNW + 2;
    LNW.iC = 2 * mNW + 3;

    mNE = sNE.b(1);
    nNE = 2 * mNE + 3;
    NNE = sum(sNE.b);
    ENE = sum(sNE.b(1:4));
    icNE = sNE.idx_corners;
    pxcNE = sNE.px(icNE);
    pycNE = sNE.py(icNE);

    jc1NE = find((abs(pxcNE - sNE.ax) < tol) & (abs(pycNE - sNE.ay) < tol));
    jc2NE = find((abs(pxcNE - sNE.bx) < tol) & (abs(pycNE - sNE.ay) < tol));
    jc3NE = find((abs(pxcNE - sNE.ax) < tol) & (abs(pycNE - sNE.by) < tol));
    jc4NE = find((abs(pxcNE - sNE.bx) < tol) & (abs(pycNE - sNE.by) < tol));

    r1NE = ENE + jc1NE;
    r2NE = ENE + jc2NE;
    r3NE = ENE + jc3NE;
    r4NE = ENE + jc4NE;

    ii0NE = zeros(NNE, 1);
    BNE = sparse(NNE, nNE);
    FNE = sparse(nNE, nNE);
    gNE = zeros(nNE, 1);

    ii0NE(mNE + (1:mNE)) = IBC{2}(sNE.px(sNE.idx_right), sNE.py(sNE.idx_right));
    ii0NE(3 * mNE + (1:mNE)) = IBC{4}(sNE.px(sNE.idx_top), sNE.py(sNE.idx_top));
    ii0NE(r2NE) = 0.5 * IBC{2}(sNE.px(icNE(jc2NE)), sNE.py(icNE(jc2NE)));
    ii0NE(r3NE) = 0.5 * IBC{4}(sNE.px(icNE(jc3NE)), sNE.py(icNE(jc3NE)));
    ii0NE(r4NE) = 0.5 * (IBC{2}(sNE.px(icNE(jc4NE)), sNE.py(icNE(jc4NE))) + IBC{4}(sNE.px(icNE(jc4NE)), sNE.py(icNE(jc4NE))));

    BNE(1:mNE, 1:mNE) = speye(mNE);
    BNE(2 * mNE + (1:mNE), mNE + (1:mNE)) = speye(mNE);
    BNE(r2NE, 2 * mNE + 1) = 0.5;
    BNE(r3NE, 2 * mNE + 2) = 0.5;
    BNE(r1NE, 2 * mNE + 3) = 1.0;

    FNE(2 * mNE + 1, 2 * mNE + 1) = 0.5;
    FNE(2 * mNE + 2, 2 * mNE + 2) = 0.5;
    gNE(2 * mNE + 1) = -0.5 * IBC{2}(sNE.px(icNE(jc2NE)), sNE.py(icNE(jc2NE)));
    gNE(2 * mNE + 2) = -0.5 * IBC{4}(sNE.px(icNE(jc3NE)), sNE.py(icNE(jc3NE)));

    rowsNE = [1:mNE, 2 * mNE + (1:mNE), r2NE, r3NE, r1NE];

    LNE.T = HNE(rowsNE, :) * BNE + FNE;
    LNE.h = HNE(rowsNE, :) * ii0NE + PNE(rowsNE) + gNE;
    LNE.B = BNE;
    LNE.ii0 = ii0NE;
    LNE.n = nNE;
    LNE.m = mNE;
    LNE.mc = full(sNE.Mb(rowsNE(end), rowsNE(end)));
    LNE.iB = mNE + (1:mNE);
    LNE.iRM = 2 * mNE + 1;
    LNE.iTM = 2 * mNE + 2;
    LNE.iC = 2 * mNE + 3;

    eSW_BR = sparse(1, nSW);
    eSW_BR(LSW.iBR) = 1;

    eSE_BM = sparse(1, nSE);
    eSE_BM(LSE.iBM) = 1;

    eSW_TL = sparse(1, nSW);
    eSW_TL(LSW.iTL) = 1;

    eNW_LM = sparse(1, nNW);
    eNW_LM(LNW.iLM) = 1;

    eNW_TM = sparse(1, nNW);
    eNW_TM(LNW.iTM) = 1;

    eNE_TM = sparse(1, nNE);
    eNE_TM(LNE.iTM) = 1;

    eSE_RM = sparse(1, nSE);
    eSE_RM(LSE.iRM) = 1;

    eNE_RM = sparse(1, nNE);
    eNE_RM(LNE.iRM) = 1;

    eSW_C = sparse(1, nSW);
    eSW_C(LSW.iC) = 1;

    eSE_C = sparse(1, nSE);
    eSE_C(LSE.iC) = 1;

    eNW_C = sparse(1, nNW);
    eNW_C(LNW.iC) = 1;

    eNE_C = sparse(1, nNE);
    eNE_C(LNE.iC) = 1;

    ZmSW = sparse(mSW, nSW);
    ZmSE = sparse(mSW, nSE);
    ZmNW = sparse(mSW, nNW);
    ZmNE = sparse(mSW, nNE);

    Z1SW = sparse(1, nSW);
    Z1SE = sparse(1, nSE);
    Z1NW = sparse(1, nNW);
    Z1NE = sparse(1, nNE);

    ISW_R = ZmSW;
    ISW_R(:, LSW.iR) = speye(mSW);

    ISE_L = ZmSE;
    ISE_L(:, 1:mSE) = speye(mSE);

    INW_R = ZmNW;
    INW_R(:, LNW.iR) = speye(mNW);

    INE_L = ZmNE;
    INE_L(:, 1:mNE) = speye(mNE);

    ISW_T = ZmSW;
    ISW_T(:, LSW.iT) = speye(mSW);

    INW_B = ZmNW;
    INW_B(:, LNW.iB) = speye(mNW);

    ISE_T = ZmSE;
    ISE_T(:, LSE.iT) = speye(mSE);

    INE_B = ZmNE;
    INE_B(:, LNE.iB) = speye(mNE);

    A = [
        ISW_R,                               LSE.T(1:mSW,:),                        ZmNW,                               ZmNE;
        LSW.T(LSW.iR,:),                     ISE_L,                                 ZmNW,                               ZmNE;
        ZmSW,                                ZmSE,                                  INW_R,                              LNE.T(1:mNW,:);
        ZmSW,                                ZmSE,                                  LNW.T(LNW.iR,:),                    INE_L;
        ISW_T,                               ZmSE,                                  LNW.T(LNW.iB,:),                    ZmNE;
        LSW.T(LSW.iT,:),                     ZmSE,                                  INW_B,                              ZmNE;
        ZmSW,                                ISE_T,                                 ZmNW,                               LNE.T(LNE.iB,:);
        ZmSW,                                LSE.T(LSE.iT,:),                       ZmNW,                               INE_B;
        eSW_BR,                              LSE.T(LSE.iBM,:),                      Z1NW,                               Z1NE;
        LSW.T(LSW.iBR,:),                    eSE_BM,                                Z1NW,                               Z1NE;
        eSW_TL,                              Z1SE,                                  LNW.T(LNW.iLM,:),                   Z1NE;
        LSW.T(LSW.iTL,:),                    Z1SE,                                  eNW_LM,                             Z1NE;
        Z1SW,                                Z1SE,                                  eNW_TM,                             LNE.T(LNE.iTM,:);
        Z1SW,                                Z1SE,                                  LNW.T(LNW.iTM,:),                   eNE_TM;
        Z1SW,                                eSE_RM,                                Z1NW,                               LNE.T(LNE.iRM,:);
        Z1SW,                                LSE.T(LSE.iRM,:),                      Z1NW,                               eNE_RM;
        (eSW_C - LSW.T(LSW.iC,:)) / LSW.mc,  (-eSE_C + LSE.T(LSE.iC,:)) / LSE.mc,   Z1NW,                               Z1NE;
        (eSW_C - LSW.T(LSW.iC,:)) / LSW.mc,  Z1SE,                                  (-eNW_C + LNW.T(LNW.iC,:)) / LNW.mc, Z1NE;
        (eSW_C - LSW.T(LSW.iC,:)) / LSW.mc,  Z1SE,                                  Z1NW,                               (-eNE_C + LNE.T(LNE.iC,:)) / LNE.mc;
        eSW_C + LSW.T(LSW.iC,:),             eSE_C + LSE.T(LSE.iC,:),               eNW_C + LNW.T(LNW.iC,:),            eNE_C + LNE.T(LNE.iC,:)
    ];

    b = [
        -LSE.h(1:mSW);
        -LSW.h(LSW.iR);
        -LNE.h(1:mNW);
        -LNW.h(LNW.iR);
        -LNW.h(LNW.iB);
        -LSW.h(LSW.iT);
        -LNE.h(LNE.iB);
        -LSE.h(LSE.iT);
        -LSE.h(LSE.iBM);
        -LSW.h(LSW.iBR);
        -LNW.h(LNW.iLM);
        -LSW.h(LSW.iTL);
        -LNE.h(LNE.iTM);
        -LNW.h(LNW.iTM);
        -LNE.h(LNE.iRM);
        -LSE.h(LSE.iRM);
        LSW.h(LSW.iC) / LSW.mc - LSE.h(LSE.iC) / LSE.mc;
        LSW.h(LSW.iC) / LSW.mc - LNW.h(LNW.iC) / LNW.mc;
        LSW.h(LSW.iC) / LSW.mc - LNE.h(LNE.iC) / LNE.mc;
        -(LSW.h(LSW.iC) + LSE.h(LSE.iC) + LNW.h(LNW.iC) + LNE.h(LNE.iC))
    ];

    x = A \ b;

    I_SW = 1:nSW;
    I_SE = nSW + (1:nSE);
    I_NW = nSW + nSE + (1:nNW);
    I_NE = nSW + nSE + nNW + (1:nNE);

    xSW = x(I_SW);
    xSE = x(I_SE);
    xNW = x(I_NW);
    xNE = x(I_NE);

    iiSW = LSW.ii0 + LSW.B * xSW;
    iiSE = LSE.ii0 + LSE.B * xSE;
    iiNW = LNW.ii0 + LNW.B * xNW;
    iiNE = LNE.ii0 + LNE.B * xNE;

    ooSW = HSW * iiSW + PSW;
    ooSE = HSE * iiSE + PSE;
    ooNW = HNW * iiNW + PNW;
    ooNE = HNE * iiNE + PNE;

    ubSW = (iiSW - ooSW) / (2i * eta);
    ubSE = (iiSE - ooSE) / (2i * eta);
    ubNW = (iiNW - ooNW) / (2i * eta);
    ubNE = (iiNE - ooNE) / (2i * eta);

    uSW = zeros(size(sSW.rhs));
    uSE = zeros(size(sSE.rhs));
    uNW = zeros(size(sNW.rhs));
    uNE = zeros(size(sNE.rhs));

    uSW(sSW.idx_interior) = sSW.A \ (sSW.rhs(sSW.idx_interior) - sSW.B * ubSW);
    uSE(sSE.idx_interior) = sSE.A \ (sSE.rhs(sSE.idx_interior) - sSE.B * ubSE);
    uNW(sNW.idx_interior) = sNW.A \ (sNW.rhs(sNW.idx_interior) - sNW.B * ubNW);
    uNE(sNE.idx_interior) = sNE.A \ (sNE.rhs(sNE.idx_interior) - sNE.B * ubNE);

    uSW(sSW.idx_boundary) = ubSW;
    uSE(sSE.idx_boundary) = ubSE;
    uNW(sNW.idx_boundary) = ubNW;
    uNE(sNE.idx_boundary) = ubNE;

    uexSW = u_exact(sSW.px, sSW.py);
    uexSE = u_exact(sSE.px, sSE.py);
    uexNW = u_exact(sNW.px, sNW.py);
    uexNE = u_exact(sNE.px, sNE.py);

    num = norm(uSW - uexSW)^2 + norm(uSE - uexSE)^2 + norm(uNW - uexNW)^2 + norm(uNE - uexNE)^2;
    den = norm(uexSW)^2 + norm(uexSE)^2 + norm(uexNW)^2 + norm(uexNE)^2;

    fprintf("div: %i\n", div);
    fprintf("Rel L2 error = %.15e\n", sqrt(num / den));

    n = 2^div + 1;

    figure;
    hold on;
    surf(reshape(sSW.px, n, n), reshape(sSW.py, n, n), reshape(real(uSW), n, n));
    surf(reshape(sSE.px, n, n), reshape(sSE.py, n, n), reshape(real(uSE), n, n));
    surf(reshape(sNW.px, n, n), reshape(sNW.py, n, n), reshape(real(uNW), n, n));
    surf(reshape(sNE.px, n, n), reshape(sNE.py, n, n), reshape(real(uNE), n, n));
    view(3);
    hold off;
end