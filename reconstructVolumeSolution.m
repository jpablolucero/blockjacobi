function [u_global] = reconstructVolumeSolution(s,divP,div,u_ref,nd,u_skel)
b  = 2^div - 1;
m  = size(nd.elementsPerFace,1) * b;
Nx = (2^divP) * (2^div) + 1;

uI = cell(numel(s),1);
for e = 1:numel(s)
    uB = zeros(sum(s{e}.b),1);

    p = 0;
    for f = 1:4
        k = nd.facePerElement(e,f);
        if k ~= 0
            uB(p+(1:b)) = u_skel((k-1)*b + (1:b));
        else
            uB(p+(1:b)) = u_ref(s{e}.px(s{e}.idx(f)), s{e}.py(s{e}.idx(f)));
        end
        p = p + b;
    end

    for j = 1:4
        kp = nd.crossPointsPerElement(e,j);
        if kp ~= 0
            uB(p+j) = u_skel(m + kp);
        else
            uB(p+j) = u_ref(s{e}.px(s{e}.idx_corners(j)), s{e}.py(s{e}.idx_corners(j)));
        end
    end

    uI{e} = s{e}.A \ (s{e}.rhs_interior - s{e}.B * uB);
end

rowsA = 2^div + 1;
u_global = zeros(Nx^2,1);

P = 2^divP;
for e = 1:numel(s)
    uB = zeros(sum(s{e}.b),1);

    p = 0;
    for f = 1:4
        k = nd.facePerElement(e,f);
        if k ~= 0
            uB(p+(1:b)) = u_skel((k-1)*b + (1:b));
        else
            uB(p+(1:b)) = u_ref(s{e}.px(s{e}.idx(f)), s{e}.py(s{e}.idx(f)));
        end
        p = p + b;
    end
    for j = 1:4
        kp = nd.crossPointsPerElement(e,j);
        if kp ~= 0
            uB(p+j) = u_skel(m + kp);
        else
            uB(p+j) = u_ref(s{e}.px(s{e}.idx_corners(j)), s{e}.py(s{e}.idx_corners(j)));
        end
    end

    u_loc = zeros(numel(s{e}.px),1);
    u_loc(s{e}.idx_interior) = s{e}.A \ (s{e}.rhs_interior - s{e}.B * uB);
    u_loc(s{e}.idx_boundary) = uB;

    [il, jl] = ind2sub([rowsA rowsA], (1:numel(s{e}.px))');
    ix = mod(e-1, P) + 1;
    jy = floor((e-1)/P) + 1;

    gx = (ix-1)*(rowsA-1) + il;
    gy = (jy-1)*(rowsA-1) + jl;
    ids = (gy-1)*Nx + gx;

    u_global(ids) = u_loc;
end
end
