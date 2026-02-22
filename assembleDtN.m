function [S,R] = assembleDtN(u_ref,s,div,nd)
dofPerFace       = 2^div-1;
totalFaceDofs    = size(nd.elementsPerFace,1)*dofPerFace;
totalCrossPoints = size(nd.crossPointGrid,1);

S = spalloc(totalFaceDofs+totalCrossPoints, ...
            totalFaceDofs+totalCrossPoints, 0); 
R = sparse(totalFaceDofs+totalCrossPoints,1);

for i = 1:size(nd.elementsPerFace,1)
    row0 = (i-1)*dofPerFace;
    R(row0+(1:dofPerFace)) = s{nd.elementsPerFace(i,2)}.h{nd.elementSidePerFace(i,2)} ...
                           + s{nd.elementsPerFace(i,1)}.h{nd.elementSidePerFace(i,1)};
    for f = 1:4

        k = nd.facePerElement(nd.elementsPerFace(i,1),f); 
        T = s{nd.elementsPerFace(i,1)}.T{nd.elementSidePerFace(i,1),f};
        if k~=0
            S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) = ...
                S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) + T;
        else
            R(row0+(1:dofPerFace)) = R(row0+(1:dofPerFace)) - ...
                                     T*u_ref(s{nd.elementsPerFace(i,1)}.px(s{nd.elementsPerFace(i,1)}.idx(f)), ...
                                             s{nd.elementsPerFace(i,1)}.py(s{nd.elementsPerFace(i,1)}.idx(f)));
        end

        k = nd.facePerElement(nd.elementsPerFace(i,2),f); 
        T = s{nd.elementsPerFace(i,2)}.T{nd.elementSidePerFace(i,2),f};
        if k~=0
            S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) = ...
                S(row0+(1:dofPerFace),(k-1)*dofPerFace+(1:dofPerFace)) + T;
        else
            R(row0+(1:dofPerFace)) = R(row0+(1:dofPerFace)) - ...
                                     T*u_ref(s{nd.elementsPerFace(i,2)}.px(s{nd.elementsPerFace(i,2)}.idx(f)), ...
                                             s{nd.elementsPerFace(i,2)}.py(s{nd.elementsPerFace(i,2)}.idx(f)));
        end
    end
    for c = 1:4
        
        k = nd.crossPointsPerElement(nd.elementsPerFace(i,1),c);
        T = s{nd.elementsPerFace(i,1)}.T{nd.elementSidePerFace(i,1),5}(:,c);
        if k~=0
            S(row0+(1:dofPerFace), totalFaceDofs+k) = ...
                S(row0+(1:dofPerFace), totalFaceDofs+k) + T;
        else
            R(row0+(1:dofPerFace)) = R(row0+(1:dofPerFace)) - ...
                                     T*u_ref(s{nd.elementsPerFace(i,1)}.px(s{nd.elementsPerFace(i,1)}.idx_corners(c)), ...
                                             s{nd.elementsPerFace(i,1)}.py(s{nd.elementsPerFace(i,1)}.idx_corners(c)));
        end
        
        k = nd.crossPointsPerElement(nd.elementsPerFace(i,2),c);
        T = s{nd.elementsPerFace(i,2)}.T{nd.elementSidePerFace(i,2),5}(:,c);
        if k~=0
            S(row0+(1:dofPerFace), totalFaceDofs+k) = ...
                S(row0+(1:dofPerFace), totalFaceDofs+k) + T;
        else
            R(row0+(1:dofPerFace)) = R(row0+(1:dofPerFace)) - ...
                                     T*u_ref(s{nd.elementsPerFace(i,2)}.px(s{nd.elementsPerFace(i,2)}.idx_corners(c)), ...
                                             s{nd.elementsPerFace(i,2)}.py(s{nd.elementsPerFace(i,2)}.idx_corners(c)));
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
                lc = ic; break;
            end
        end

        for f = 1:4
            kk = nd.facePerElement(e,f);
            T5f = s{e}.T{5,f}(lc,:);
            if kk ~= 0
                col0 = (kk-1)*dofPerFace;
                S(r, col0+(1:dofPerFace)) = S(r, col0+(1:dofPerFace)) + T5f;
            else
                R(r) = R(r) - T5f * u_ref( ...
                    s{e}.px(s{e}.idx(f)), s{e}.py(s{e}.idx(f)) );
            end
        end

        for jc = 1:4
            kp = nd.crossPointsPerElement(e,jc);
            a = s{e}.T{5,5}(lc,jc);
            if kp ~= 0
                S(r, totalFaceDofs+kp) = S(r, totalFaceDofs+kp) + a;
            else
                R(r) = R(r) - a * u_ref( ...
                    s{e}.px(s{e}.idx_corners(jc)), ...
                    s{e}.py(s{e}.idx_corners(jc)) );
            end
        end

        R(r) = R(r) + s{e}.h{5}(lc);
    end
end

end
