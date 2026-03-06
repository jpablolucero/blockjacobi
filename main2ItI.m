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
    % Subdomains. W: West, E: East.
    sW = Subdomain(div, rhs, axW, bxW, ay, by, 0, c0, "DtN", @get_sem);
    sE = Subdomain(div, rhs, axE, bxE, ay, by, 0, c0, "DtN", @get_sem);

    % Solution values at the corner of each subdomain
    uCornerW_zg = BC{5}(sW.px(sW.idx_corners), sW.py(sW.idx_corners));
    uCornerE_zg = BC{5}(sE.px(sE.idx_corners), sE.py(sE.idx_corners));

    % Solution values at global boundary set to zero, interface corners
    % kept, zg : Zeroed-out Global corners
    uCornerW_zg([1,3]) = 0;
    uCornerE_zg([2,4]) = 0;

    % Amount of DoFs at edges without corners for each subdomain
    allEdgesSizeW = sum(sW.b(1:4));
    allEdgesSizeE = sum(sE.b(1:4));

    % Homogeneous Dirichlet-to-Neumann operator for each subdomain
    DtNHW = cell2mat(sW.T(1:5,1:5));
    DtNHE = cell2mat(sE.T(1:5,1:5));

    % Index sets of boundaries 
    % wic: Without Interface Corners
    idx_wic_W = [(1:allEdgesSizeW), allEdgesSizeW+1, allEdgesSizeW+3].';
    idx_wic_E = [(1:allEdgesSizeE), allEdgesSizeE+2, allEdgesSizeE+4].';

    % Homogeneous Dirichlet-to-Neumann operator 
    % wic: Without Interface Corners
    DtNH_wic_W = DtNHW(idx_wic_W,idx_wic_W);
    DtNH_wic_E = DtNHE(idx_wic_E,idx_wic_E);

    % Particular Dirichlet-to-Neumann operator with
    % zg : Zeroed-out Global corners
    DtNP_zg_W = cell2mat(sW.h(1:5)) - DtNHW*[zeros(allEdgesSizeW,1); uCornerW_zg];
    DtNP_zg_E = cell2mat(sE.h(1:5)) - DtNHE*[zeros(allEdgesSizeE,1); uCornerE_zg];

    % Particular Dirichlet-to-Neumann operator with 
    % zg : Zeroed-out Global corners,
    % wic: Without Interface Corners
    DtNP_wic_zg_W = DtNP_zg_W(idx_wic_W);
    DtNP_wic_zg_E = DtNP_zg_E(idx_wic_E);

    % Homogeneous Impedance-to-Impedance operator with 
    % zg : Zeroed-out Global corners,
    % wic: Without Interface Corners
    ItIH_wic_zg_W = sW.Mb(idx_wic_W,idx_wic_W) \ (DtNH_wic_W - 1i*eta*sW.Mb(idx_wic_W,idx_wic_W)) * ((DtNH_wic_W + 1i*eta*sW.Mb(idx_wic_W,idx_wic_W)) \ sW.Mb(idx_wic_W,idx_wic_W)) ;
    ItIH_wic_zg_E = sE.Mb(idx_wic_E,idx_wic_E) \ (DtNH_wic_E - 1i*eta*sE.Mb(idx_wic_E,idx_wic_E)) * ((DtNH_wic_E + 1i*eta*sE.Mb(idx_wic_E,idx_wic_E)) \ sE.Mb(idx_wic_E,idx_wic_E)) ;

    % Particular Impedance-to-Impedance operator with
    % zg : Zeroed-out Global corners,
    % wic: Without Interface Corners
    ItIP_wic_zg_W = -(2i*eta) * ((DtNH_wic_W + 1i*eta*sW.Mb(idx_wic_W,idx_wic_W)) \ DtNP_wic_zg_W);
    ItIP_wic_zg_E = -(2i*eta) * ((DtNH_wic_E + 1i*eta*sE.Mb(idx_wic_E,idx_wic_E)) \ DtNP_wic_zg_E);

    % Cell Homogeneous Impedance-to-Impedance operator with
    % zg : Zeroed-out Global corners,
    % wic: Without Interface Corners
    ItIH_wic_zg_W_cell = mat2cell(ItIH_wic_zg_W, [sW.b(1:4), 2], [sW.b(1:4), 2]);  
    ItIH_wic_zg_E_cell = mat2cell(ItIP_wic_zg_W, [sW.b(1:4), 2], 1);

    % Cell Particular Impedance-to-Impedance operator with
    % zg : Zeroed-out Global corners,
    % wic: Without Interface Corners
    ItIP_wic_zg_W_cell = mat2cell(ItIH_wic_zg_E, [sE.b(1:4), 2], [sE.b(1:4), 2]);  
    ItIP_wic_zg_E_cell = mat2cell(ItIP_wic_zg_E, [sE.b(1:4), 2], 1);
    
    % West subdomain non-interface corner points
    xW1 = sW.px(sW.idx_corners(1));  yW1 = sW.py(sW.idx_corners(1));
    xW3 = sW.px(sW.idx_corners(3));  yW3 = sW.py(sW.idx_corners(3));

    % Incoming impedance boundary condition 
    iiW = {IBC{1}(sW.px(sW.idx_left),    sW.py(sW.idx_left)),...
           zeros(sW.b(2),1),...
           IBC{3}(sW.px(sW.idx_bottom),  sW.py(sW.idx_bottom)),...
           IBC{4}(sW.px(sW.idx_top),     sW.py(sW.idx_top)),...
           [0.5*(IBC{1}(xW1,yW1) + IBC{3}(xW1,yW1));...
            0.5*(IBC{1}(xW3,yW3) + IBC{4}(xW3,yW3))];};

    xE2 = sE.px(sE.idx_corners(2));  yE2 = sE.py(sE.idx_corners(2));
    xE4 = sE.px(sE.idx_corners(4));  yE4 = sE.py(sE.idx_corners(4));
    iiE = {zeros(sE.b(1),1),...
           IBC{2}(sE.px(sE.idx_right),   sE.py(sE.idx_right)),...
           IBC{3}(sE.px(sE.idx_bottom),  sE.py(sE.idx_bottom)),...
           IBC{4}(sE.px(sE.idx_top),     sE.py(sE.idx_top)),...
           [0.5*(IBC{2}(xE2,yE2) + IBC{3}(xE2,yE2));...
            0.5*(IBC{2}(xE4,yE4) + IBC{4}(xE4,yE4))]};

    rW = ItIH_wic_zg_W_cell{2,1}*iiW{1} + ItIH_wic_zg_W_cell{2,3}*iiW{3} + ItIH_wic_zg_W_cell{2,4}*iiW{4} + ItIH_wic_zg_W_cell{2,5}*iiW{5} + ItIH_wic_zg_E_cell{2};
    rE = ItIP_wic_zg_W_cell{1,2}*iiE{2} + ItIP_wic_zg_W_cell{1,3}*iiE{3} + ItIP_wic_zg_W_cell{1,4}*iiE{4} + ItIP_wic_zg_W_cell{1,5}*iiE{5} + ItIP_wic_zg_E_cell{1};

    m  = sW.b(2);
    oI = [speye(m), ItIH_wic_zg_W_cell{2,2};
          ItIP_wic_zg_W_cell{1,1},  speye(m)] \ [rW; rE];

    iiW{2} = -oI(m+(1:m));
    iiE{1} = -oI(1:m);

    ooW = cell(5,1);
    ooE = cell(5,1);

    ooW{1} = ItIH_wic_zg_W_cell{1,1}*iiW{1} + ItIH_wic_zg_W_cell{1,2}*iiW{2} + ItIH_wic_zg_W_cell{1,3}*iiW{3} + ItIH_wic_zg_W_cell{1,4}*iiW{4} + ItIH_wic_zg_W_cell{1,5}*iiW{5} + ItIH_wic_zg_E_cell{1};
    ooW{2} = ItIH_wic_zg_W_cell{2,1}*iiW{1} + ItIH_wic_zg_W_cell{2,2}*iiW{2} + ItIH_wic_zg_W_cell{2,3}*iiW{3} + ItIH_wic_zg_W_cell{2,4}*iiW{4} + ItIH_wic_zg_W_cell{2,5}*iiW{5} + ItIH_wic_zg_E_cell{2};
    ooW{3} = ItIH_wic_zg_W_cell{3,1}*iiW{1} + ItIH_wic_zg_W_cell{3,2}*iiW{2} + ItIH_wic_zg_W_cell{3,3}*iiW{3} + ItIH_wic_zg_W_cell{3,4}*iiW{4} + ItIH_wic_zg_W_cell{3,5}*iiW{5} + ItIH_wic_zg_E_cell{3};
    ooW{4} = ItIH_wic_zg_W_cell{4,1}*iiW{1} + ItIH_wic_zg_W_cell{4,2}*iiW{2} + ItIH_wic_zg_W_cell{4,3}*iiW{3} + ItIH_wic_zg_W_cell{4,4}*iiW{4} + ItIH_wic_zg_W_cell{4,5}*iiW{5} + ItIH_wic_zg_E_cell{4};
    ooW{5} = ItIH_wic_zg_W_cell{5,1}*iiW{1} + ItIH_wic_zg_W_cell{5,2}*iiW{2} + ItIH_wic_zg_W_cell{5,3}*iiW{3} + ItIH_wic_zg_W_cell{5,4}*iiW{4} + ItIH_wic_zg_W_cell{5,5}*iiW{5} + ItIH_wic_zg_E_cell{5};

    ooE{1} = ItIP_wic_zg_W_cell{1,1}*iiE{1} + ItIP_wic_zg_W_cell{1,2}*iiE{2} + ItIP_wic_zg_W_cell{1,3}*iiE{3} + ItIP_wic_zg_W_cell{1,4}*iiE{4} + ItIP_wic_zg_W_cell{1,5}*iiE{5} + ItIP_wic_zg_E_cell{1};
    ooE{2} = ItIP_wic_zg_W_cell{2,1}*iiE{1} + ItIP_wic_zg_W_cell{2,2}*iiE{2} + ItIP_wic_zg_W_cell{2,3}*iiE{3} + ItIP_wic_zg_W_cell{2,4}*iiE{4} + ItIP_wic_zg_W_cell{2,5}*iiE{5} + ItIP_wic_zg_E_cell{2};
    ooE{3} = ItIP_wic_zg_W_cell{3,1}*iiE{1} + ItIP_wic_zg_W_cell{3,2}*iiE{2} + ItIP_wic_zg_W_cell{3,3}*iiE{3} + ItIP_wic_zg_W_cell{3,4}*iiE{4} + ItIP_wic_zg_W_cell{3,5}*iiE{5} + ItIP_wic_zg_E_cell{3};
    ooE{4} = ItIP_wic_zg_W_cell{4,1}*iiE{1} + ItIP_wic_zg_W_cell{4,2}*iiE{2} + ItIP_wic_zg_W_cell{4,3}*iiE{3} + ItIP_wic_zg_W_cell{4,4}*iiE{4} + ItIP_wic_zg_W_cell{4,5}*iiE{5} + ItIP_wic_zg_E_cell{4};
    ooE{5} = ItIP_wic_zg_W_cell{5,1}*iiE{1} + ItIP_wic_zg_W_cell{5,2}*iiE{2} + ItIP_wic_zg_W_cell{5,3}*iiE{3} + ItIP_wic_zg_W_cell{5,4}*iiE{4} + ItIP_wic_zg_W_cell{5,5}*iiE{5} + ItIP_wic_zg_E_cell{5};

    ubW{1} = (iiW{1} - ooW{1})/(2i*eta);
    ubW{2} = (iiW{2} - ooW{2})/(2i*eta);
    ubW{3} = (iiW{3} - ooW{3})/(2i*eta);
    ubW{4} = (iiW{4} - ooW{4})/(2i*eta);

    ubE{1} = (iiE{1} - ooE{1})/(2i*eta);
    ubE{2} = (iiE{2} - ooE{2})/(2i*eta);
    ubE{3} = (iiE{3} - ooE{3})/(2i*eta);
    ubE{4} = (iiE{4} - ooE{4})/(2i*eta);

    ubW5 = BC{5}(sW.px(sW.idx_corners), sW.py(sW.idx_corners));
    ubE5 = BC{5}(sE.px(sE.idx_corners), sE.py(sE.idx_corners));

    ubW5([1,3]) = (iiW{5} - ooW{5})/(2i*eta);
    ubE5([2,4]) = (iiE{5} - ooE{5})/(2i*eta);

    ubW{5} = ubW5;
    ubE{5} = ubE5;

    ubWv = [ubW{1}; ubW{2}; ubW{3}; ubW{4}; ubW{5}];
    ubEv = [ubE{1}; ubE{2}; ubE{3}; ubE{4}; ubE{5}];

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

    figure;
    hold on;
    surf(reshape(sW.px,nW,nW), reshape(sW.py,nW,nW), reshape(real(uW),nW,nW));
    surf(reshape(sE.px,nE,nE), reshape(sE.py,nE,nE), reshape(real(uE),nE,nE));
    view(3);
    hold off;
end