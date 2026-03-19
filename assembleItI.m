function [S,R] = assembleItI(s,div,nd,IBC)
dofPerFace = 2^div - 1;
nFaces = size(nd.elementsPerFace,1);
totalFaceDofs = 2 * nFaces * dofPerFace;

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

faceBlockPerElementSide = zeros(size(nd.facePerElement));
for i = 1:nFaces
    e = nd.elementsPerFace(i,1);
    f = nd.elementSidePerFace(i,1);
    faceBlockPerElementSide(e,f) = 2*(i-1) + 1;

    e = nd.elementsPerFace(i,2);
    f = nd.elementSidePerFace(i,2);
    faceBlockPerElementSide(e,f) = 2*(i-1) + 2;
end

S = spalloc(totalFaceDofs + totalPointDofs, totalFaceDofs + totalPointDofs, 0);
R = sparse(totalFaceDofs + totalPointDofs, 1);

for i = 1:nFaces
    for q = 1:2
        block = 2*(i-1) + q;
        rows = (block-1)*dofPerFace + (1:dofPerFace);

        e = nd.elementsPerFace(i,q);
        side = nd.elementSidePerFace(i,q);
        otherBlock = 2*(i-1) + (3-q);

        R(rows) = -s{e}.h{side};

        for f = 1:4
            k = nd.facePerElement(e,f);
            if k ~= 0
                b = faceBlockPerElementSide(e,f);
                cols = (b-1)*dofPerFace + (1:dofPerFace);
                S(rows,cols) = S(rows,cols) + s{e}.T{side,f};
            end
        end

        colsOther = (otherBlock-1)*dofPerFace + (1:dofPerFace);
        S(rows,colsOther) = S(rows,colsOther) + speye(dofPerFace);

        for c = 1:4
            k = nd.crossPointsPerElement(e,c);
            if k ~= 0
                pd = pointDofPerElement(e,c);
                if nnz(nd.elementsPerCrossPoint(k,:)) == 2
                    w = 0.5;
                else
                    w = 1.0;
                end
                S(rows,totalFaceDofs + pd) = S(rows,totalFaceDofs + pd) + w * s{e}.T{side,5}(:,c);
            end
        end
    end
end

cornerMbIndex = @(lc) 4 * dofPerFace + lc;

Xmin = min(cellfun(@(t) t.ax, s));
Xmax = max(cellfun(@(t) t.bx, s));
Ymin = min(cellfun(@(t) t.ay, s));
Ymax = max(cellfun(@(t) t.by, s));
tol = 1.0e-12;
cornerExtSides = {[1,3], [2,3], [1,4], [2,4]};

for k = 1:size(nd.crossPointGrid,1)
    E = nd.elementsPerCrossPoint(k,:);
    E = E(E ~= 0);
    nInc = numel(E);

    lcs = zeros(1,nInc);
    pds = zeros(1,nInc);

    for t = 1:nInc
        e = E(t);
        for c = 1:4
            if nd.crossPointsPerElement(e,c) == k
                lcs(t) = c;
                pds(t) = pointDofPerElement(e,c);
                break
            end
        end
    end

    if nInc == 2
        for t = 1:2
            e = E(t);
            lc = lcs(t);
            pd = pds(t);
            otherPd = pds(3-t);

            r = totalFaceDofs + pd;

            for f = 1:4
                kk = nd.facePerElement(e,f);
                if kk ~= 0
                    b = faceBlockPerElementSide(e,f);
                    cols = (b-1)*dofPerFace + (1:dofPerFace);
                    S(r,cols) = S(r,cols) + s{e}.T{5,f}(lc,:);
                end
            end

            for jc = 1:4
                kp = nd.crossPointsPerElement(e,jc);
                if kp ~= 0
                    qd = pointDofPerElement(e,jc);
                    if nnz(nd.elementsPerCrossPoint(kp,:)) == 2
                        w = 0.5;
                    else
                        w = 1.0;
                    end
                    S(r,totalFaceDofs + qd) = S(r,totalFaceDofs + qd) + w * s{e}.T{5,5}(lc,jc);
                end
            end

            S(r,totalFaceDofs + pd) = S(r,totalFaceDofs + pd) + 0.5;
            S(r,totalFaceDofs + otherPd) = S(r,totalFaceDofs + otherPd) + 1.0;

            xc = s{e}.px(s{e}.idx_corners(lc));
            yc = s{e}.py(s{e}.idx_corners(lc));
            halfIBC = 0;
            for ss = cornerExtSides{lc}
                sideIsExt = (ss == 1 && abs(s{e}.ax - Xmin) < tol) || ...
                            (ss == 2 && abs(s{e}.bx - Xmax) < tol) || ...
                            (ss == 3 && abs(s{e}.ay - Ymin) < tol) || ...
                            (ss == 4 && abs(s{e}.by - Ymax) < tol);
                if sideIsExt
                    halfIBC = halfIBC + 0.5 * IBC{ss}(xc,yc);
                end
            end

            R(r) = -(s{e}.h{5}(lc) - halfIBC);
        end
    elseif nInc == 4
        eRef = E(1);
        lcRef = lcs(1);
        pdRef = pds(1);
        MRef = s{eRef}.Mb(cornerMbIndex(lcRef),cornerMbIndex(lcRef));

        for t = 2:4
            e = E(t);
            lc = lcs(t);
            pd = pds(t);
            Mt = s{e}.Mb(cornerMbIndex(lc),cornerMbIndex(lc));

            r = totalFaceDofs + pds(t-1);

            for f = 1:4
                kk = nd.facePerElement(eRef,f);
                if kk ~= 0
                    b = faceBlockPerElementSide(eRef,f);
                    cols = (b-1)*dofPerFace + (1:dofPerFace);
                    S(r,cols) = S(r,cols) - s{eRef}.T{5,f}(lcRef,:) / MRef;
                end
            end

            for jc = 1:4
                kp = nd.crossPointsPerElement(eRef,jc);
                if kp ~= 0
                    qd = pointDofPerElement(eRef,jc);
                    if nnz(nd.elementsPerCrossPoint(kp,:)) == 2
                        w = 0.5;
                    else
                        w = 1.0;
                    end
                    S(r,totalFaceDofs + qd) = S(r,totalFaceDofs + qd) - w * s{eRef}.T{5,5}(lcRef,jc) / MRef;
                end
            end

            S(r,totalFaceDofs + pdRef) = S(r,totalFaceDofs + pdRef) + 1.0 / MRef;

            for f = 1:4
                kk = nd.facePerElement(e,f);
                if kk ~= 0
                    b = faceBlockPerElementSide(e,f);
                    cols = (b-1)*dofPerFace + (1:dofPerFace);
                    S(r,cols) = S(r,cols) + s{e}.T{5,f}(lc,:) / Mt;
                end
            end

            for jc = 1:4
                kp = nd.crossPointsPerElement(e,jc);
                if kp ~= 0
                    qd = pointDofPerElement(e,jc);
                    if nnz(nd.elementsPerCrossPoint(kp,:)) == 2
                        w = 0.5;
                    else
                        w = 1.0;
                    end
                    S(r,totalFaceDofs + qd) = S(r,totalFaceDofs + qd) + w * s{e}.T{5,5}(lc,jc) / Mt;
                end
            end

            S(r,totalFaceDofs + pd) = S(r,totalFaceDofs + pd) - 1.0 / Mt;
            R(r) = s{eRef}.h{5}(lcRef) / MRef - s{e}.h{5}(lc) / Mt;
        end

        r = totalFaceDofs + pds(4);

        for t = 1:4
            e = E(t);
            lc = lcs(t);
            pd = pds(t);

            for f = 1:4
                kk = nd.facePerElement(e,f);
                if kk ~= 0
                    b = faceBlockPerElementSide(e,f);
                    cols = (b-1)*dofPerFace + (1:dofPerFace);
                    S(r,cols) = S(r,cols) + s{e}.T{5,f}(lc,:);
                end
            end

            for jc = 1:4
                kp = nd.crossPointsPerElement(e,jc);
                if kp ~= 0
                    qd = pointDofPerElement(e,jc);
                    if nnz(nd.elementsPerCrossPoint(kp,:)) == 2
                        w = 0.5;
                    else
                        w = 1.0;
                    end
                    S(r,totalFaceDofs + qd) = S(r,totalFaceDofs + qd) + w * s{e}.T{5,5}(lc,jc);
                end
            end

            S(r,totalFaceDofs + pd) = S(r,totalFaceDofs + pd) + 1.0;
            R(r) = R(r) - s{e}.h{5}(lc);
        end
    end
end
end