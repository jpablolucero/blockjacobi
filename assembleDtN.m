function [S,R] = assembleDtN(BC,s,div,nd)
dofPerFace       = 2^div-1;
totalFaceDofs    = size(nd.elementsPerFace,1)*dofPerFace;
totalCrossPoints = size(nd.crossPointGrid,1);

S = spalloc(totalFaceDofs+totalCrossPoints, totalFaceDofs+totalCrossPoints, 0);
R = sparse(totalFaceDofs+totalCrossPoints,1);

tol = 1.0e-14;

for i = 1:size(nd.elementsPerFace,1)
    row0 = (i-1)*dofPerFace;

    R(row0+(1:dofPerFace)) = s{nd.elementsPerFace(i,2)}.h{nd.elementSidePerFace(i,2)} ...
                           + s{nd.elementsPerFace(i,1)}.h{nd.elementSidePerFace(i,1)};

    for f = 1:4
        e = nd.elementsPerFace(i,1);
        k = nd.facePerElement(e,f);
        T = s{e}.T{nd.elementSidePerFace(i,1),f};
        if k~=0
            S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) = ...
                S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) + T;
        else
            x = s{e}.px(s{e}.idx(f));
            y = s{e}.py(s{e}.idx(f));
            R(row0+(1:dofPerFace)) = R(row0+(1:dofPerFace)) - T*BC{f}(x,y);
        end

        e = nd.elementsPerFace(i,2);
        k = nd.facePerElement(e,f);
        T = s{e}.T{nd.elementSidePerFace(i,2),f};
        if k~=0
            S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) = ...
                S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) + T;
        else
            x = s{e}.px(s{e}.idx(f));
            y = s{e}.py(s{e}.idx(f));
            R(row0+(1:dofPerFace)) = R(row0+(1:dofPerFace)) - T*BC{f}(x,y);
        end
    end

    for c = 1:4
        e = nd.elementsPerFace(i,1);
        k = nd.crossPointsPerElement(e,c);
        T = s{e}.T{nd.elementSidePerFace(i,1),5}(:,c);
        if k~=0
            S(row0+(1:dofPerFace), totalFaceDofs+k) = ...
                S(row0+(1:dofPerFace), totalFaceDofs+k) + T;
        else
            xc = s{e}.px(s{e}.idx_corners(c));
            yc = s{e}.py(s{e}.idx_corners(c));
            uc = evalCornerBC(BC, xc, yc, s{e}.px, s{e}.py, tol);
            R(row0+(1:dofPerFace)) = R(row0+(1:dofPerFace)) - T*uc;
        end

        e = nd.elementsPerFace(i,2);
        k = nd.crossPointsPerElement(e,c);
        T = s{e}.T{nd.elementSidePerFace(i,2),5}(:,c);
        if k~=0
            S(row0+(1:dofPerFace), totalFaceDofs+k) = ...
                S(row0+(1:dofPerFace), totalFaceDofs+k) + T;
        else
            xc = s{e}.px(s{e}.idx_corners(c));
            yc = s{e}.py(s{e}.idx_corners(c));
            uc = evalCornerBC(BC, xc, yc, s{e}.px, s{e}.py, tol);
            R(row0+(1:dofPerFace)) = R(row0+(1:dofPerFace)) - T*uc;
        end
    end
end

for k = 1:totalCrossPoints
    r = totalFaceDofs + k;
    E = nd.elementsPerCrossPoint(k,:);

    for t = 1:4
        e = E(t);

        lc = 0;
        for ic = 1:4
            if nd.crossPointsPerElement(e,ic) == k
                lc = ic;
                break;
            end
        end

        for f = 1:4
            kk = nd.facePerElement(e,f);
            T5f = s{e}.T{5,f}(lc,:);
            if kk ~= 0
                col0 = (kk-1)*dofPerFace;
                S(r, col0+(1:dofPerFace)) = S(r, col0+(1:dofPerFace)) + T5f;
            else
                x = s{e}.px(s{e}.idx(f));
                y = s{e}.py(s{e}.idx(f));
                R(r) = R(r) - T5f * BC{f}(x,y);
            end
        end

        for jc = 1:4
            kp = nd.crossPointsPerElement(e,jc);
            a = s{e}.T{5,5}(lc,jc);
            if kp ~= 0
                S(r, totalFaceDofs+kp) = S(r, totalFaceDofs+kp) + a;
            else
                xc = s{e}.px(s{e}.idx_corners(jc));
                yc = s{e}.py(s{e}.idx_corners(jc));
                uc = evalCornerBC(BC, xc, yc, s{e}.px, s{e}.py, tol);
                R(r) = R(r) - a * uc;
            end
        end

        R(r) = R(r) + s{e}.h{5}(lc);
    end
end

    function u = evalCornerBC(BC, x, y, px, py, tol)
        ax = min(px);
        bx = max(px);
        ay = min(py);
        by = max(py);

        on_left   = abs(x-ax) < tol;
        on_right  = abs(x-bx) < tol;
        on_bottom = abs(y-ay) < tol;
        on_top    = abs(y-by) < tol;

        vals = [];

        if on_left
            vals(end+1,1) = BC{1}(x,y);
        end
        if on_right
            vals(end+1,1) = BC{2}(x,y);
        end
        if on_bottom
            vals(end+1,1) = BC{3}(x,y);
        end
        if on_top
            vals(end+1,1) = BC{4}(x,y);
        end

        if isempty(vals)
            error('assembleDtN: Corner point not detected as lying on any boundary.');
        end

        u = mean(vals);
    end

end