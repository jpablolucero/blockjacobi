function [u_global, u_cells, pxG, pyG] = reconstructVolumeSolutionItI( ...
    s, divP, div, u_ref, u_skel, skel, IBC, eta)
%RECONSTRUCTVOLUMESOLUTIONITI  Recover the volume solution from skeleton data.
%
%   [u_global, u_cells, pxG, pyG] = reconstructVolumeSolutionItI(
%       s, divP, div, u_ref, u_skel, skel, IBC, eta)
%
%   Inputs
%     s       – cell array of Subdomain objects (linear index, jy-major)
%     divP    – partition level
%     div     – element refinement level
%     u_ref   – exact / reference solution handle  u_ref(px,py)
%     u_skel  – skeleton solution vector  (column)
%     skel    – struct returned by assembleItI (numbering data)
%     IBC     – cell{4} of impedance boundary-condition functions
%     eta     – impedance parameter
%
%   Outputs
%     u_global – global solution vector  (column, nG*nG x 1)
%     u_cells  – cell array of per-subdomain solution vectors
%     pxG      – global x-coordinates   (column, nG*nG x 1)
%     pyG      – global y-coordinates   (column, nG*nG x 1)

nSub        = skel.nSub;
m           = skel.m;
edgeStart   = skel.edgeStart;
cornerGIdx  = skel.cornerGIdx;
cornerWeight = skel.cornerWeight;
cornerBCval = skel.cornerBCval;
isIntMat    = skel.isIntMat;

Sub = reshape(s, [nSub, nSub]);

n  = 2^div + 1;                          % nodes per subdomain side
nG = nSub * (n - 1) + 1;                 % global nodes per direction

u_cells = cell(nSub, nSub);
pxG     = zeros(nG * nG, 1);
pyG     = zeros(nG * nG, 1);
u_global = zeros(nG * nG, 1);

for jj = 1:nSub
    for ii = 1:nSub
        Si = Sub{ii,jj};

        % ---- Build the incoming-trace vector  (4m + 4 entries) ----
        incoming = zeros(4*m + 4, 1);
        for ss = 1:4
            idx = (ss-1)*m + (1:m);
            if isIntMat(ii,jj,ss)
                incoming(idx) = u_skel(edgeStart(ii,jj,ss) : edgeStart(ii,jj,ss)+m-1);
            end
        end
        for c = 1:4
            w = cornerWeight(ii,jj,c);
            if w == 0, continue; end
            incoming(4*m + c) = w * u_skel(cornerGIdx(ii,jj,c));
        end

        % ---- Outgoing = T * incoming + h ----
        oo = cell2mat(Si.T) * incoming + cell2mat(Si.h);

        % ---- Build the full boundary trace ----
        fullTrace = zeros(4*m + 4, 1);
        for ss = 1:4
            idx = (ss-1)*m + (1:m);
            if isIntMat(ii,jj,ss)
                fullTrace(idx) = u_skel(edgeStart(ii,jj,ss) : edgeStart(ii,jj,ss)+m-1);
            else
                switch ss
                    case 1, fullTrace(idx) = IBC{1}(Si.px(Si.idx_left),   Si.py(Si.idx_left));
                    case 2, fullTrace(idx) = IBC{2}(Si.px(Si.idx_right),  Si.py(Si.idx_right));
                    case 3, fullTrace(idx) = IBC{3}(Si.px(Si.idx_bottom), Si.py(Si.idx_bottom));
                    case 4, fullTrace(idx) = IBC{4}(Si.px(Si.idx_top),    Si.py(Si.idx_top));
                end
            end
        end
        for c = 1:4
            w   = cornerWeight(ii,jj,c);
            bcv = cornerBCval{ii,jj}(c);
            if w == 0
                fullTrace(4*m + c) = bcv;
            else
                fullTrace(4*m + c) = bcv + w * u_skel(cornerGIdx(ii,jj,c));
            end
        end

        % ---- Recover physical boundary values ----
        ub = (fullTrace - oo) / (2i * eta);

        % ---- Interior solve ----
        u = zeros(size(Si.rhs));
        u(Si.idx_interior) = Si.A \ (Si.rhs(Si.idx_interior) - Si.B * ub);
        u(Si.idx_boundary) = ub;

        u_cells{ii,jj} = u;

        % ---- Place into global grid (vectorised) ----
        %  Local node (lx, ly) in 1..n maps to global (gx, gy):
        %    gx = (ii-1)*(n-1) + lx,   gy = (jj-1)*(n-1) + ly
        %  Both local and global arrays are column-major.
        lx = (1:n).';                                     % local x-indices
        ly = (1:n).';                                     % local y-indices
        gx = (ii-1)*(n-1) + lx;                           % global x-indices
        gy = (jj-1)*(n-1) + ly;                           % global y-indices

        locIdx = reshape((ly'-1)*n + lx, [], 1);          % n*n local indices
        gloIdx = reshape((gy'-1)*nG + gx, [], 1);         % n*n global indices

        pxG(gloIdx)      = Si.px(locIdx);
        pyG(gloIdx)      = Si.py(locIdx);
        u_global(gloIdx) = u(locIdx);
    end
end

end