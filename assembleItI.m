function [S, R, skel] = assembleItI(s, div, divP, IBC, Xdom, Ydom)
%ASSEMBLEITI  Assemble the Impedance-to-Impedance skeleton system.
%
%   [S, R, skel] = assembleItI(s, div, divP, IBC, Xdom, Ydom)
%
%   Inputs
%     s      – cell array of Subdomain objects (linear index, jy-major)
%     div    – element refinement level  (2^div+1 nodes per subdomain side)
%     divP   – partition level           (2^divP subdomains per direction)
%     IBC    – cell{4} of impedance boundary-condition functions
%     Xdom   – [x0, x1] domain x-extent
%     Ydom   – [y0, y1] domain y-extent
%
%   Outputs
%     S, R   – skeleton linear system  S * u_skel = R
%     skel   – struct with numbering data needed by reconstructVolumeSolutionItI

nSub = 2^divP;
m    = 2^div - 1;

cBL = 4*m+1;  cBR = 4*m+2;  cTL = 4*m+3;  cTR = 4*m+4;
localCornerBdryIdx = [cBL, cBR, cTL, cTR];

% Reshape linear cell array to (ix, jy) grid
Sub = reshape(s, [nSub, nSub]);

% ================================================================
%  Geometry tables
% ================================================================
%  Neighbour of side s:  side nbSide(s) of subdomain offset by (nbDi,nbDj)
nbDi   = [-1,  1,  0,  0];
nbDj   = [ 0,  0, -1,  1];
nbSide = [ 2,  1,  4,  3];

%  Corner c of (ii,jj) sits at vertex (ii+cDi(c), jj+cDj(c))
cDi = [-1,  0, -1,  0];   % BL BR TL TR
cDj = [-1, -1,  0,  0];

%  Which physical-boundary sides can contribute an IBC at corner c?
cornerExtSides = {[1,3], [2,3], [1,4], [2,4]};

%  Is side s of subdomain (ii,jj) an interior interface?
isIntMat = false(nSub, nSub, 4);
for jj = 1:nSub
    for ii = 1:nSub
        for ss = 1:4
            isIntMat(ii,jj,ss) = ~( (ss==1 && ii==1)    || (ss==2 && ii==nSub) || ...
                                    (ss==3 && jj==1)     || (ss==4 && jj==nSub) );
        end
    end
end

% ================================================================
%  Corner weights and BC values
% ================================================================
cornerWeight = zeros(nSub, nSub, 4);
cornerBCval  = cell(nSub, nSub);

for jj = 1:nSub
    for ii = 1:nSub
        bcv = zeros(4,1);
        for c = 1:4
            vi = ii + cDi(c);
            vj = jj + cDj(c);
            nExtEdge = (vi==0) + (vi==nSub) + (vj==0) + (vj==nSub);

            if     nExtEdge >= 2,  cornerWeight(ii,jj,c) = 0;
            elseif nExtEdge == 1,  cornerWeight(ii,jj,c) = 0.5;
            else,                  cornerWeight(ii,jj,c) = 1;
            end

            cidx = localCornerBdryIdx(c);
            px_c = Sub{ii,jj}.px(Sub{ii,jj}.idx_boundary(cidx));
            py_c = Sub{ii,jj}.py(Sub{ii,jj}.idx_boundary(cidx));
            for ss = cornerExtSides{c}(:).'
                if ~isIntMat(ii,jj,ss)
                    bcv(c) = bcv(c) + 0.5 * IBC{ss}(px_c, py_c);
                end
            end
        end
        cornerBCval{ii,jj} = bcv;
    end
end

% ================================================================
%  Number the skeleton unknowns
% ================================================================
%  --- Edge unknowns: m per interior side ---
edgeStart = zeros(nSub, nSub, 4);
cnt = 0;
for jj = 1:nSub
    for ii = 1:nSub
        for ss = 1:4
            if isIntMat(ii,jj,ss)
                edgeStart(ii,jj,ss) = cnt*m + 1;
                cnt = cnt + 1;
            end
        end
    end
end
nEdgeUnk = cnt * m;

%  --- Corner unknowns: one per subdomain-corner at each non-domain-corner vertex ---
cornerGIdx = zeros(nSub, nSub, 4);
nCornerUnk = 0;

%  Vertex (vi,vj) ∈ {0,...,nSub}^2 is shared by up to 4 subdomains:
%    TR(c=4) of (vi,  vj  )     if vi>=1,    vj>=1
%    TL(c=3) of (vi+1,vj  )     if vi<nSub,  vj>=1
%    BR(c=2) of (vi,  vj+1)     if vi>=1,    vj<nSub
%    BL(c=1) of (vi+1,vj+1)     if vi<nSub,  vj<nSub
vertexList = struct('vi',{},'vj',{},'subs',{},'type',{});

for vj = 0:nSub
    for vi = 0:nSub
        nExtEdge = (vi==0) + (vi==nSub) + (vj==0) + (vj==nSub);
        if nExtEdge >= 2, continue; end      % domain corner – fully determined

        subs = zeros(0,3);
        if vi>=1      && vj>=1,      subs(end+1,:) = [vi,   vj,   4]; end %#ok<AGROW>
        if vi+1<=nSub && vj>=1,      subs(end+1,:) = [vi+1, vj,   3]; end %#ok<AGROW>
        if vi>=1      && vj+1<=nSub, subs(end+1,:) = [vi,   vj+1, 2]; end %#ok<AGROW>
        if vi+1<=nSub && vj+1<=nSub, subs(end+1,:) = [vi+1, vj+1, 1]; end %#ok<AGROW>

        for kk = 1:size(subs,1)
            nCornerUnk = nCornerUnk + 1;
            cornerGIdx(subs(kk,1), subs(kk,2), subs(kk,3)) = nEdgeUnk + nCornerUnk;
        end

        if nExtEdge == 1, vtype = 'junction';
        else,             vtype = 'interior';
        end
        vertexList(end+1) = struct('vi',vi,'vj',vj,'subs',subs,'type',vtype); %#ok<AGROW>
    end
end
N = nEdgeUnk + nCornerUnk;

% ================================================================
%  Assemble skeleton system  S x = R
% ================================================================
estNNZ = nEdgeUnk * (2*m + 10) + nCornerUnk * (2*m + 16);
II = zeros(estNNZ,1);  JJ = zeros(estNNZ,1);  VV = zeros(estNNZ,1);
ptr = 0;
R   = zeros(N,1);

% Helper to grow triplet arrays if needed
    function ensureCapacity(needed)
        if ptr + needed > length(II)
            extra = max(estNNZ, needed);
            II = [II; zeros(extra,1)]; %#ok<SETNU>
            JJ = [JJ; zeros(extra,1)]; %#ok<SETNU>
            VV = [VV; zeros(extra,1)]; %#ok<SETNU>
        end
    end

% -------- Edge equations --------
for jj = 1:nSub
    for ii = 1:nSub
        Si = Sub{ii,jj};
        for ss = 1:4
            if ~isIntMat(ii,jj,ss), continue; end
            rowOff = edgeStart(ii,jj,ss);
            rows   = (rowOff : rowOff+m-1).';

            % T{s,t}  for every interior side t of this subdomain
            for t = 1:4
                if ~isIntMat(ii,jj,t), continue; end
                colOff = edgeStart(ii,jj,t);
                cols   = (colOff : colOff+m-1).';
                [ri, ci, vi] = find(sparse(Si.T{ss,t}));
                nz = length(vi);
                ensureCapacity(nz);
                II(ptr+1:ptr+nz) = rows(ri);
                JJ(ptr+1:ptr+nz) = cols(ci);
                VV(ptr+1:ptr+nz) = vi;
                ptr = ptr + nz;
            end

            % Identity coupling to the neighbour's matching side
            ni = ii + nbDi(ss);  nj = jj + nbDj(ss);  ns = nbSide(ss);
            colOff = edgeStart(ni, nj, ns);
            cols   = (colOff : colOff+m-1).';
            ensureCapacity(m);
            II(ptr+1:ptr+m) = rows;
            JJ(ptr+1:ptr+m) = cols;
            VV(ptr+1:ptr+m) = ones(m,1);
            ptr = ptr + m;

            % Corner coupling:  weight(c) * T{s,5}(:,c)
            for c = 1:4
                w = cornerWeight(ii,jj,c);
                if w == 0, continue; end
                gidx = cornerGIdx(ii,jj,c);
                if gidx == 0, continue; end
                vals = full(w * Si.T{ss,5}(:,c));
                idx  = find(vals);  nz = length(idx);
                ensureCapacity(nz);
                II(ptr+1:ptr+nz) = rows(idx);
                JJ(ptr+1:ptr+nz) = gidx;
                VV(ptr+1:ptr+nz) = vals(idx);
                ptr = ptr + nz;
            end

            R(rows) = -Si.h{ss};
        end
    end
end

% -------- Corner equations --------
for vv = 1:length(vertexList)
    V    = vertexList(vv);
    subs = V.subs;
    nv   = size(subs,1);

    if strcmp(V.type, 'junction')
        % ==== Junction vertex (2 subdomains, 2 equations) ====
        for kk = 1:nv
            ii = subs(kk,1);  jj = subs(kk,2);  c = subs(kk,3);
            Si  = Sub{ii,jj};
            row = cornerGIdx(ii,jj,c);

            % Edge terms: T{5,t}(c,:)
            for t = 1:4
                if ~isIntMat(ii,jj,t), continue; end
                colOff = edgeStart(ii,jj,t);
                vals   = full(Si.T{5,t}(c,:));
                idx    = find(vals);  nz = length(idx);
                ensureCapacity(nz);
                II(ptr+1:ptr+nz) = row;
                JJ(ptr+1:ptr+nz) = colOff - 1 + idx(:);
                VV(ptr+1:ptr+nz) = vals(idx(:));
                ptr = ptr + nz;
            end

            % Corner terms:  w*T{5,5}(c,c2) + delta(c,c2)*w
            for c2 = 1:4
                w2   = cornerWeight(ii,jj,c2);
                if w2 == 0, continue; end
                gidx = cornerGIdx(ii,jj,c2);
                if gidx == 0, continue; end
                val  = w2 * Si.T{5,5}(c,c2);
                if c2 == c, val = val + w2; end
                if val ~= 0
                    ensureCapacity(1);
                    ptr = ptr+1;
                    II(ptr) = row; JJ(ptr) = gidx; VV(ptr) = val;
                end
            end

            % Identity coupling (=1) to each other subdomain's corner at this vertex
            for kk2 = 1:nv
                if kk2 == kk, continue; end
                gidx = cornerGIdx(subs(kk2,1), subs(kk2,2), subs(kk2,3));
                ensureCapacity(1);
                ptr = ptr+1;
                II(ptr) = row; JJ(ptr) = gidx; VV(ptr) = 1;
            end

            R(row) = -(Si.h{5}(c) - cornerBCval{ii,jj}(c));
        end

    else
        % ==== Interior vertex (4 subdomains, 4 equations) ====
        %  Reference = first entry.
        %    ref row   ->  sum equation
        %    other rows ->  difference equations (ref - other) / Mb
        ii_r = subs(1,1); jj_r = subs(1,2); c_r = subs(1,3);
        S_r  = Sub{ii_r, jj_r};
        Mb_r = S_r.Mb(localCornerBdryIdx(c_r), localCornerBdryIdx(c_r));

        % --- Difference equations ---
        for other = 2:nv
            ii_o = subs(other,1); jj_o = subs(other,2); c_o = subs(other,3);
            S_o  = Sub{ii_o, jj_o};
            Mb_o = S_o.Mb(localCornerBdryIdx(c_o), localCornerBdryIdx(c_o));
            row  = cornerGIdx(ii_o, jj_o, c_o);

            % Edges from reference  (-T{5,t}(c_r,:) / Mb_r)
            for t = 1:4
                if ~isIntMat(ii_r,jj_r,t), continue; end
                colOff = edgeStart(ii_r,jj_r,t);
                vals   = full(-S_r.T{5,t}(c_r,:) / Mb_r);
                idx    = find(vals);  nz = length(idx);
                ensureCapacity(nz);
                II(ptr+1:ptr+nz) = row;
                JJ(ptr+1:ptr+nz) = colOff - 1 + idx(:);
                VV(ptr+1:ptr+nz) = vals(idx(:));
                ptr = ptr + nz;
            end

            % Edges from other  (+T{5,t}(c_o,:) / Mb_o)
            for t = 1:4
                if ~isIntMat(ii_o,jj_o,t), continue; end
                colOff = edgeStart(ii_o,jj_o,t);
                vals   = full(S_o.T{5,t}(c_o,:) / Mb_o);
                idx    = find(vals);  nz = length(idx);
                ensureCapacity(nz);
                II(ptr+1:ptr+nz) = row;
                JJ(ptr+1:ptr+nz) = colOff - 1 + idx(:);
                VV(ptr+1:ptr+nz) = vals(idx(:));
                ptr = ptr + nz;
            end

            % Corners from reference: -w*T{5,5}(c_r,c2)/Mb_r  (+w/Mb_r if c2==c_r)
            for c2 = 1:4
                w2   = cornerWeight(ii_r,jj_r,c2);
                if w2 == 0, continue; end
                gidx = cornerGIdx(ii_r,jj_r,c2);
                if gidx == 0, continue; end
                val  = -w2 * S_r.T{5,5}(c_r,c2) / Mb_r;
                if c2 == c_r, val = val + w2 / Mb_r; end
                if val ~= 0
                    ensureCapacity(1);
                    ptr = ptr+1;
                    II(ptr) = row; JJ(ptr) = gidx; VV(ptr) = val;
                end
            end

            % Corners from other: +w*T{5,5}(c_o,c2)/Mb_o  (-w/Mb_o if c2==c_o)
            for c2 = 1:4
                w2   = cornerWeight(ii_o,jj_o,c2);
                if w2 == 0, continue; end
                gidx = cornerGIdx(ii_o,jj_o,c2);
                if gidx == 0, continue; end
                val  = w2 * S_o.T{5,5}(c_o,c2) / Mb_o;
                if c2 == c_o, val = val - w2 / Mb_o; end
                if val ~= 0
                    ensureCapacity(1);
                    ptr = ptr+1;
                    II(ptr) = row; JJ(ptr) = gidx; VV(ptr) = val;
                end
            end

            R(row) = S_r.h{5}(c_r)/Mb_r - S_o.h{5}(c_o)/Mb_o;
        end

        % --- Sum equation (placed in the reference corner's row) ---
        row = cornerGIdx(ii_r, jj_r, c_r);

        for kk = 1:nv
            ii_k = subs(kk,1); jj_k = subs(kk,2); c_k = subs(kk,3);
            S_k  = Sub{ii_k, jj_k};

            % Edges: +T{5,t}(c_k,:)
            for t = 1:4
                if ~isIntMat(ii_k,jj_k,t), continue; end
                colOff = edgeStart(ii_k,jj_k,t);
                vals   = full(S_k.T{5,t}(c_k,:));
                idx    = find(vals);  nz = length(idx);
                ensureCapacity(nz);
                II(ptr+1:ptr+nz) = row;
                JJ(ptr+1:ptr+nz) = colOff - 1 + idx(:);
                VV(ptr+1:ptr+nz) = vals(idx(:));
                ptr = ptr + nz;
            end

            % Corners: w*T{5,5}(c_k,c2) + delta*w
            for c2 = 1:4
                w2   = cornerWeight(ii_k,jj_k,c2);
                if w2 == 0, continue; end
                gidx = cornerGIdx(ii_k,jj_k,c2);
                if gidx == 0, continue; end
                val  = w2 * S_k.T{5,5}(c_k,c2);
                if c2 == c_k, val = val + w2; end
                if val ~= 0
                    ensureCapacity(1);
                    ptr = ptr+1;
                    II(ptr) = row; JJ(ptr) = gidx; VV(ptr) = val;
                end
            end
        end

        hsum = 0;
        for kk = 1:nv
            hsum = hsum + Sub{subs(kk,1),subs(kk,2)}.h{5}(subs(kk,3));
        end
        R(row) = -hsum;
    end
end

% Build sparse matrix
II = II(1:ptr);  JJ = JJ(1:ptr);  VV = VV(1:ptr);
S  = sparse(II, JJ, VV, N, N);

% Pack skeleton numbering data for reconstruction
skel.edgeStart          = edgeStart;
skel.cornerGIdx         = cornerGIdx;
skel.cornerWeight       = cornerWeight;
skel.cornerBCval        = cornerBCval;
skel.isIntMat           = isIntMat;
skel.localCornerBdryIdx = localCornerBdryIdx;
skel.m                  = m;
skel.nSub               = nSub;

end