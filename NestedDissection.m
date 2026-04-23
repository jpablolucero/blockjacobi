classdef NestedDissection < handle

    properties
        N
        map
        elementsPerFace
        elementSidePerFace
        facePerElement
        levels
        crossPointGrid
        crossPointsPerElement
        elementsPerCrossPoint
        enclosedCrossPointsPerMacroFace
        macroFacesPerLevel
        permutation
        nDofsPerMacroFace
        separatorBlockSizes
        separatorPermutation
    end

    methods
        function obj = NestedDissection(arg)
            if nargin == 0
                return
            end

            if isscalar(arg)
                obj.N = arg;
                n = 2^obj.N;
                A = reshape(1:(n * n), [n, n]).';
            else
                A = arg;
                obj.N = log2(size(arg, 1));
            end

            obj.map = cell(1, 1);
            obj.map{1, 1} = cell(1, 1);
            obj.map{1, 1}{1, 1} = A;

            obj.elementsPerFace = zeros(0, 2);
            obj.elementSidePerFace = zeros(0, 2);
            obj.facePerElement = zeros(2^obj.N * 2^obj.N, 4);

            obj.levels = [];

            obj.crossPointGrid = [];
            for j = 1:2^obj.N - 1
                for i = 1:2^obj.N - 1
                    obj.crossPointGrid(end + 1, :) = [i, j];
                end
            end

            obj.crossPointsPerElement = zeros(2^obj.N * 2^obj.N, 4);
            obj.elementsPerCrossPoint = zeros((2^obj.N - 1) * (2^obj.N - 1), 4);

            for i = 1:(2^obj.N - 1) * (2^obj.N - 1)
                obj.crossPointsPerElement(A(obj.crossPointGrid(i, 2), obj.crossPointGrid(i, 1)), 4) = i;
                obj.crossPointsPerElement(A(obj.crossPointGrid(i, 2), obj.crossPointGrid(i, 1) + 1), 3) = i;
                obj.crossPointsPerElement(A(obj.crossPointGrid(i, 2) + 1, obj.crossPointGrid(i, 1)), 2) = i;
                obj.crossPointsPerElement(A(obj.crossPointGrid(i, 2) + 1, obj.crossPointGrid(i, 1) + 1), 1) = i;

                obj.elementsPerCrossPoint(i, 1) = A(obj.crossPointGrid(i, 2), obj.crossPointGrid(i, 1));
                obj.elementsPerCrossPoint(i, 2) = A(obj.crossPointGrid(i, 2), obj.crossPointGrid(i, 1) + 1);
                obj.elementsPerCrossPoint(i, 3) = A(obj.crossPointGrid(i, 2) + 1, obj.crossPointGrid(i, 1));
                obj.elementsPerCrossPoint(i, 4) = A(obj.crossPointGrid(i, 2) + 1, obj.crossPointGrid(i, 1) + 1);
            end

            obj.enclosedCrossPointsPerMacroFace = {};
            obj.macroFacesPerLevel = zeros(1, 0);

            obj.permutation = zeros(0, 1);
            obj.nDofsPerMacroFace = {};
            obj.separatorBlockSizes = {};
            obj.separatorPermutation = {};
        end

        function divide(obj, flag)

            prevGrid = obj.map{end, 1};
            A0 = prevGrid{1, 1};

            if flag == 0
                if size(A0, 1) < 2
                    return
                end
                nextGrid = cell(2 * size(prevGrid, 1), size(prevGrid, 2));
            else
                if size(A0, 2) < 2
                    return
                end
                nextGrid = cell(size(prevGrid, 1), 2 * size(prevGrid, 2));
            end

            facePairs = zeros(0, 2);
            mergeFace = zeros(0, 2);

            if flag == 0
                verticesPerMergedFace = zeros(numel(prevGrid), size(A0, 2) - 1);
            else
                verticesPerMergedFace = zeros(numel(prevGrid), size(A0, 1) - 1);
            end

            if flag == 0
                for c = 1:size(prevGrid, 2)
                    for r = 1:size(prevGrid, 1)
                        A = prevGrid{r, c};
                        mid = size(A, 1) / 2;
                        nextGrid{2 * r - 1, c} = A(1:mid, :);
                        nextGrid{2 * r, c} = A(mid + 1:end, :);

                        for j = 1:size(A, 2)
                            facePairs(end + 1, :) = [A(mid, j), A(mid + 1, j)];
                            mergeFace(end + 1, :) = [4, 3];
                        end

                        for j = 1:size(A, 2) - 1
                            verticesPerMergedFace((c - 1) * size(prevGrid, 1) + r, j) = obj.crossPointsPerElement(A(mid, j), 4);
                        end
                    end
                end
            else
                for r = 1:size(prevGrid, 1)
                    for c = 1:size(prevGrid, 2)
                        A = prevGrid{r, c};
                        mid = size(A, 2) / 2;
                        nextGrid{r, 2 * c - 1} = A(:, 1:mid);
                        nextGrid{r, 2 * c} = A(:, mid + 1:end);

                        for i = 1:size(A, 1)
                            facePairs(end + 1, :) = [A(i, mid), A(i, mid + 1)];
                            mergeFace(end + 1, :) = [2, 1];
                        end

                        for i = 1:size(A, 1) - 1
                            verticesPerMergedFace((r - 1) * size(prevGrid, 2) + c, i) = obj.crossPointsPerElement(A(i, mid), 4);
                        end
                    end
                end
            end

            obj.map{end + 1, 1} = nextGrid;

            obj.divide(1 - flag);

            for i = 1:size(facePairs, 1)
                if flag == 1
                    obj.facePerElement(facePairs(i, 1), 2) = size(obj.elementsPerFace, 1) + i;
                    obj.facePerElement(facePairs(i, 2), 1) = size(obj.elementsPerFace, 1) + i;
                else
                    obj.facePerElement(facePairs(i, 1), 4) = size(obj.elementsPerFace, 1) + i;
                    obj.facePerElement(facePairs(i, 2), 3) = size(obj.elementsPerFace, 1) + i;
                end
            end

            obj.levels(end + 1) = size(facePairs, 1);
            obj.elementsPerFace = [obj.elementsPerFace; facePairs];
            obj.elementSidePerFace = [obj.elementSidePerFace; mergeFace];

            obj.enclosedCrossPointsPerMacroFace{1, end + 1} = verticesPerMergedFace;
            obj.macroFacesPerLevel(end + 1) = size(verticesPerMergedFace, 1);
        end

        function calculateReorderingCross(obj, dofsPerFace)

            A = obj.map{1, 1}{1, 1};

            n = size(A, 1);

            nFaces = size(obj.elementsPerFace, 1);
            nFaceDofsTotal = nFaces * dofsPerFace;
            nCrossPoints = size(obj.crossPointGrid, 1);

            cpIndex = zeros(n - 1, n - 1);
            for cp = 1:nCrossPoints
                cpIndex(obj.crossPointGrid(cp, 2), obj.crossPointGrid(cp, 1)) = cp;
            end

            obj.permutation = zeros(0, 1);
            obj.nDofsPerMacroFace = cell(obj.N, 1);

            for level = 1:obj.N

                branchLength = 2^(level - 1);

                centersX = branchLength:(2 * branchLength):(n - 1);
                centersY = branchLength:(2 * branchLength):(n - 1);

                obj.nDofsPerMacroFace{level} = zeros(numel(centersX) * numel(centersY), 1);

                block = 0;

                for cy = centersY
                    for cx = centersX

                        block = block + 1;
                        blockStart = numel(obj.permutation) + 1;

                        for branch = 1:4

                            if branch == 1
                                stackType = 1;
                                stackStart = cy - branchLength + 1;
                                stackLen = branchLength;
                                stackFixed = cx;
                                stackCp = 0;
                            elseif branch == 2
                                stackType = 1;
                                stackStart = cy + 1;
                                stackLen = branchLength;
                                stackFixed = cx;
                                stackCp = 0;
                            elseif branch == 3
                                stackType = 2;
                                stackStart = cx - branchLength + 1;
                                stackLen = branchLength;
                                stackFixed = cy;
                                stackCp = 0;
                            else
                                stackType = 2;
                                stackStart = cx + 1;
                                stackLen = branchLength;
                                stackFixed = cy;
                                stackCp = 0;
                            end

                            while ~isempty(stackType)

                                type = stackType(end);
                                start = stackStart(end);
                                len = stackLen(end);
                                fixed = stackFixed(end);
                                cp = stackCp(end);

                                stackType(end) = [];
                                stackStart(end) = [];
                                stackLen(end) = [];
                                stackFixed(end) = [];
                                stackCp(end) = [];

                                if type == 3

                                elseif len == 1

                                    if type == 1
                                        f = obj.facePerElement(A(start, fixed), 2);
                                    else
                                        f = obj.facePerElement(A(fixed, start), 4);
                                    end

                                    faceDofs = (f - 1) * dofsPerFace + (1:dofsPerFace);
                                    obj.permutation = [obj.permutation; faceDofs(:)];

                                else

                                    half = len / 2;

                                    if type == 1
                                        cp = cpIndex(start + half - 1, fixed);
                                    else
                                        cp = cpIndex(fixed, start + half - 1);
                                    end

                                    stackType(end + 1) = type;
                                    stackStart(end + 1) = start + half;
                                    stackLen(end + 1) = half;
                                    stackFixed(end + 1) = fixed;
                                    stackCp(end + 1) = 0;

                                    stackType(end + 1) = 3;
                                    stackStart(end + 1) = 0;
                                    stackLen(end + 1) = 0;
                                    stackFixed(end + 1) = 0;
                                    stackCp(end + 1) = cp;

                                    stackType(end + 1) = type;
                                    stackStart(end + 1) = start;
                                    stackLen(end + 1) = half;
                                    stackFixed(end + 1) = fixed;
                                    stackCp(end + 1) = 0;

                                end

                            end

                        end

                        obj.nDofsPerMacroFace{level}(block) = numel(obj.permutation) - blockStart + 1;

                    end
                end
            end

            for cy = 1:(n - 1)
                for cx = 1:(n - 1)
                    obj.permutation = [obj.permutation; nFaceDofsTotal + cpIndex(cy, cx)];
                end
            end

            obj.separatorBlockSizes = obj.nDofsPerMacroFace;
            obj.separatorPermutation = cell(obj.N, 1);

            if obj.N == 1
                obj.separatorPermutation{1} = obj.permutation;
                return
            end

            gPrev = obj.getReorderD(1);
            obj.separatorPermutation{1} = obj.permutation(gPrev);

            for level = 2:obj.N

                if level == obj.N
                    gCurr = (1:numel(obj.permutation)).';
                else
                    gCurr = obj.getReorderD(level);
                end

                prefixLength = 0;
                for j = 1:(level - 1)
                    prefixLength = prefixLength + sum(obj.separatorBlockSizes{j});
                end

                prevTail = gPrev(prefixLength + 1:end);
                currTail = gCurr(prefixLength + 1:end);

                inversePrevTail = zeros(numel(obj.permutation), 1);
                for i = 1:numel(prevTail)
                    inversePrevTail(prevTail(i)) = i;
                end

                p = zeros(numel(currTail), 1);
                for i = 1:numel(currTail)
                    p(i) = inversePrevTail(currTail(i));
                end

                obj.separatorPermutation{level} = p;
                gPrev = gCurr;
            end

        end

        function permutation = getReorderD(obj, level)

            if nargin < 2
                level = 1;
            end

            A = obj.map{1, 1}{1, 1};

            n = size(A, 1);

            nFaces = size(obj.elementsPerFace, 1);
            nCrossPoints = size(obj.crossPointGrid, 1);
            nTotal = numel(obj.permutation);
            nFaceDofsTotal = nTotal - nCrossPoints;
            dofsPerFace = nFaceDofsTotal / nFaces;

            if level >= obj.N
                permutation = (1:nTotal).';
                return
            end

            branchLength = 2^(level - 1);

            prefixLength = 0;
            for j = 1:level
                prefixLength = prefixLength + sum(obj.nDofsPerMacroFace{j});
            end

            inversePermutation = zeros(nTotal, 1);
            for i = 1:nTotal
                inversePermutation(obj.permutation(i)) = i;
            end

            protectedFaces = false(nFaces, 1);

            for lev = 1:level

                h = 2^(lev - 1);

                centersX = h:(2 * h):(n - 1);
                centersY = h:(2 * h):(n - 1);

                for cy = centersY
                    for cx = centersX

                        for r = (cy - h + 1):cy
                            protectedFaces(obj.facePerElement(A(r, cx), 2)) = true;
                        end

                        for r = (cy + 1):(cy + h)
                            protectedFaces(obj.facePerElement(A(r, cx), 2)) = true;
                        end

                        for c = (cx - h + 1):cx
                            protectedFaces(obj.facePerElement(A(cy, c), 4)) = true;
                        end

                        for c = (cx + 1):(cx + h)
                            protectedFaces(obj.facePerElement(A(cy, c), 4)) = true;
                        end

                    end
                end

            end

            permutation = (1:prefixLength).';

            usedDFaces = false(nFaces, 1);

            for cy = 0:(2 * branchLength):n
                for cx = 0:(2 * branchLength):n

                    if (cx == 0 || cx == n) && (cy == 0 || cy == n)
                        continue
                    end

                    if cx > 0 && cx < n && cy > 0 && cy < n

                        faces = zeros(4 * branchLength, 1);
                        cnt = 0;

                        for r = (cy - branchLength + 1):cy
                            cnt = cnt + 1;
                            faces(cnt) = obj.facePerElement(A(r, cx), 2);
                        end

                        for r = (cy + 1):(cy + branchLength)
                            cnt = cnt + 1;
                            faces(cnt) = obj.facePerElement(A(r, cx), 2);
                        end

                        for c = (cx - branchLength + 1):cx
                            cnt = cnt + 1;
                            faces(cnt) = obj.facePerElement(A(cy, c), 4);
                        end

                        for c = (cx + 1):(cx + branchLength)
                            cnt = cnt + 1;
                            faces(cnt) = obj.facePerElement(A(cy, c), 4);
                        end

                    elseif cy == 0

                        faces = zeros(branchLength, 1);
                        cnt = 0;

                        for r = 1:branchLength
                            cnt = cnt + 1;
                            faces(cnt) = obj.facePerElement(A(r, cx), 2);
                        end

                    elseif cy == n

                        faces = zeros(branchLength, 1);
                        cnt = 0;

                        for r = (n - branchLength + 1):n
                            cnt = cnt + 1;
                            faces(cnt) = obj.facePerElement(A(r, cx), 2);
                        end

                    elseif cx == 0

                        faces = zeros(branchLength, 1);
                        cnt = 0;

                        for c = 1:branchLength
                            cnt = cnt + 1;
                            faces(cnt) = obj.facePerElement(A(cy, c), 4);
                        end

                    else

                        faces = zeros(branchLength, 1);
                        cnt = 0;

                        for c = (n - branchLength + 1):n
                            cnt = cnt + 1;
                            faces(cnt) = obj.facePerElement(A(cy, c), 4);
                        end

                    end

                    for j = 1:numel(faces)

                        if protectedFaces(faces(j))
                            error('getReorderD selected a face that belongs to an A-cross already kept in place');
                        end

                        if usedDFaces(faces(j))
                            error('getReorderD selected the same D-face twice');
                        end

                        usedDFaces(faces(j)) = true;

                        faceDofs = (faces(j) - 1) * dofsPerFace + (1:dofsPerFace);
                        for k = 1:dofsPerFace
                            permutation(end + 1, 1) = inversePermutation(faceDofs(k));
                        end

                    end

                end
            end

            if any(~protectedFaces & ~usedDFaces)
                error('getReorderD did not cover all non-protected faces');
            end

            for i = (nFaceDofsTotal + 1):nTotal
                permutation(end + 1, 1) = inversePermutation(i);
            end

        end

        function calculateReordering(obj, dofsPerFace)

            nFaceDofsTotal = size(obj.elementsPerFace, 1) * dofsPerFace;

            obj.permutation = zeros(0, 1);
            obj.nDofsPerMacroFace = cell(numel(obj.levels), 1);
            for level = 1:numel(obj.levels)
                obj.nDofsPerMacroFace{level} = zeros(obj.macroFacesPerLevel(level), 1);
            end

            faceDofCursor = 0;
            for level = 1:numel(obj.levels)
                facesPerMacroFace = obj.levels(level) / obj.macroFacesPerLevel(level);

                for macroFace = 1:obj.macroFacesPerLevel(level)
                    blockStart = numel(obj.permutation) + 1;

                    cps = obj.enclosedCrossPointsPerMacroFace{level}(macroFace, :);
                    nCps = numel(cps);

                    for f = 1:facesPerMacroFace
                        faceDofs = (faceDofCursor + (f - 1) * dofsPerFace + 1):(faceDofCursor + f * dofsPerFace);
                        obj.permutation = [obj.permutation; faceDofs(:)];

                        if f <= nCps
                            crossPointDof = nFaceDofsTotal + cps(f);
                            obj.permutation = [obj.permutation; crossPointDof];
                        end
                    end

                    obj.nDofsPerMacroFace{level}(macroFace, 1) = numel(obj.permutation) - blockStart + 1;
                    faceDofCursor = faceDofCursor + facesPerMacroFace * dofsPerFace;
                end
            end

            obj.separatorBlockSizes = obj.nDofsPerMacroFace;
            obj.separatorPermutation = cell(numel(obj.levels), 1);

            obj.separatorPermutation{1} = obj.permutation;

            nRemaining = numel(obj.permutation) - sum(obj.separatorBlockSizes{1});
            for level = 2:numel(obj.levels)
                obj.separatorPermutation{level} = (1:nRemaining).';
                nRemaining = nRemaining - sum(obj.separatorBlockSizes{level});
            end
        end

        function permutation = getSeparatorPermutation(obj, level)
            permutation = obj.separatorPermutation{level};
        end

    end
end