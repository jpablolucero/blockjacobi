function [u_global, u_cells, px_global, py_global] = reconstructVolumeSolutionItI(s,divP,div,IBC,nd,u_skel)
b = 2^div - 1;
nFaces = size(nd.elementsPerFace,1);
totalFaceDofs = 2 * nFaces * b;
Nx = (2^divP) * (2^div) + 1;

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

faceBlockPerElementSide = zeros(size(nd.facePerElement));
for i = 1:nFaces
    e = nd.elementsPerFace(i,1);
    f = nd.elementSidePerFace(i,1);
    faceBlockPerElementSide(e,f) = 2*(i-1) + 1;

    e = nd.elementsPerFace(i,2);
    f = nd.elementSidePerFace(i,2);
    faceBlockPerElementSide(e,f) = 2*(i-1) + 2;
end

rowsA = 2^div + 1;
u_global = zeros(Nx^2,1);
u_cells = cell(numel(s),1);
px_global = nan(Nx^2,1);
py_global = nan(Nx^2,1);

P = 2^divP;
Xmin = min(cellfun(@(t) t.ax, s));
Xmax = max(cellfun(@(t) t.bx, s));
Ymin = min(cellfun(@(t) t.ay, s));
Ymax = max(cellfun(@(t) t.by, s));
tol = 1.0e-12;
cornerExtSides = {[1,3], [2,3], [1,4], [2,4]};

for e = 1:numel(s)
    uInOp = zeros(sum(s{e}.b),1);
    uInFull = zeros(sum(s{e}.b),1);

    p = 0;
    for f = 1:4
        k = nd.facePerElement(e,f);
        if k ~= 0
            block = faceBlockPerElementSide(e,f);
            cols = (block-1)*b + (1:b);
            val = u_skel(cols);
            uInOp(p + (1:b)) = val;
            uInFull(p + (1:b)) = val;
        else
            uInFull(p + (1:b)) = IBC{f}(s{e}.px(s{e}.idx(f)), s{e}.py(s{e}.idx(f)));
        end
        p = p + b;
    end

    for c = 1:4
        valOp = 0;
        valFull = 0;

        xc = s{e}.px(s{e}.idx_corners(c));
        yc = s{e}.py(s{e}.idx_corners(c));

        for ss = cornerExtSides{c}
            sideIsExt = (ss == 1 && abs(s{e}.ax - Xmin) < tol) || ...
                        (ss == 2 && abs(s{e}.bx - Xmax) < tol) || ...
                        (ss == 3 && abs(s{e}.ay - Ymin) < tol) || ...
                        (ss == 4 && abs(s{e}.by - Ymax) < tol);
            if sideIsExt
                valFull = valFull + 0.5 * IBC{ss}(xc,yc);
            end
        end

        kp = nd.crossPointsPerElement(e,c);
        if kp ~= 0
            pd = pointDofPerElement(e,c);
            if nnz(nd.elementsPerCrossPoint(kp,:)) == 2
                valOp = valOp + 0.5 * u_skel(totalFaceDofs + pd);
                valFull = valFull + 0.5 * u_skel(totalFaceDofs + pd);
            else
                valOp = valOp + u_skel(totalFaceDofs + pd);
                valFull = valFull + u_skel(totalFaceDofs + pd);
            end
        end

        uInOp(p + c) = valOp;
        uInFull(p + c) = valFull;
    end

    uOut = cell2mat(s{e}.T) * uInOp + cell2mat(s{e}.h);
    uB = (uInFull - uOut) / (2i * s{e}.eta);

    u_loc = zeros(numel(s{e}.px),1);
    u_loc(s{e}.idx_interior) = s{e}.A \ (s{e}.rhs(s{e}.idx_interior) - s{e}.B * uB);
    u_loc(s{e}.idx_boundary) = uB;

    u_cells{e} = u_loc;

    [il, jl] = ind2sub([rowsA rowsA], (1:numel(s{e}.px))');
    ix = mod(e-1, P) + 1;
    jy = floor((e-1)/P) + 1;

    gx = (ix-1)*(rowsA-1) + il;
    gy = (jy-1)*(rowsA-1) + jl;
    ids = (gy-1)*Nx + gx;

    u_global(ids) = u_loc;
    px_global(ids) = s{e}.px(:);
    py_global(ids) = s{e}.py(:);
end
end