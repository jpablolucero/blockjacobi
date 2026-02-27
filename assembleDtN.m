function [S,R] = assembleDtN(s,div,nd)
dofPerFace       = 2^div-1;
totalFaceDofs    = size(nd.elementsPerFace,1)*dofPerFace;
totalCrossPoints = size(nd.crossPointGrid,1);

S = spalloc(totalFaceDofs+totalCrossPoints, totalFaceDofs+totalCrossPoints, 0);
R = sparse(totalFaceDofs+totalCrossPoints,1);

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
        end

        e = nd.elementsPerFace(i,2);
        k = nd.facePerElement(e,f);
        T = s{e}.T{nd.elementSidePerFace(i,2),f};
        if k~=0
            S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) = ...
                S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) + T;
        end
    end

    for c = 1:4
        e = nd.elementsPerFace(i,1);
        k = nd.crossPointsPerElement(e,c);
        T = s{e}.T{nd.elementSidePerFace(i,1),5}(:,c);
        if k~=0
            S(row0+(1:dofPerFace), totalFaceDofs+k) = ...
                S(row0+(1:dofPerFace), totalFaceDofs+k) + T;
        end

        e = nd.elementsPerFace(i,2);
        k = nd.crossPointsPerElement(e,c);
        T = s{e}.T{nd.elementSidePerFace(i,2),5}(:,c);
        if k~=0
            S(row0+(1:dofPerFace), totalFaceDofs+k) = ...
                S(row0+(1:dofPerFace), totalFaceDofs+k) + T;
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
            end
        end

        for jc = 1:4
            kp = nd.crossPointsPerElement(e,jc);
            a = s{e}.T{5,5}(lc,jc);
            if kp ~= 0
                S(r, totalFaceDofs+kp) = S(r, totalFaceDofs+kp) + a;
            end
        end

        R(r) = R(r) + s{e}.h{5}(lc);
    end
end
end