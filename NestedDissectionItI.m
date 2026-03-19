classdef NestedDissectionItI < handle

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
    end

    methods
        function obj = NestedDissectionItI(arg)
            if nargin == 0
                return
            end

            if isscalar(arg)
                obj.N = arg;
                n = 2^obj.N;
                A = reshape(1:(n*n), [n, n]).';
            else
                A = arg;
                obj.N = log2(size(arg, 1));
            end

            n = 2^obj.N;

            obj.map = cell(1,1);
            obj.map{1,1} = cell(1,1);
            obj.map{1,1}{1,1} = A;

            obj.elementsPerFace = zeros(0,2);
            obj.elementSidePerFace = zeros(0,2);
            obj.facePerElement = zeros(n*n,4);
            obj.levels = [];

            obj.crossPointGrid = [];
            for j = 0:n
                for i = 0:n
                    if ~(i == 0 && j == 0) && ~(i == 0 && j == n) && ~(i == n && j == 0) && ~(i == n && j == n)
                        obj.crossPointGrid(end+1,:) = [i,j];
                    end
                end
            end

            nCrossPoints = size(obj.crossPointGrid,1);
            obj.crossPointsPerElement = zeros(n*n,4);
            obj.elementsPerCrossPoint = zeros(nCrossPoints,4);

            for k = 1:nCrossPoints
                i = obj.crossPointGrid(k,1);
                j = obj.crossPointGrid(k,2);

                if i >= 1 && i <= n-1 && j >= 1 && j <= n-1
                    obj.crossPointsPerElement(A(j,i),4)     = k;
                    obj.crossPointsPerElement(A(j,i+1),3)   = k;
                    obj.crossPointsPerElement(A(j+1,i),2)   = k;
                    obj.crossPointsPerElement(A(j+1,i+1),1) = k;

                    obj.elementsPerCrossPoint(k,1) = A(j,i);
                    obj.elementsPerCrossPoint(k,2) = A(j,i+1);
                    obj.elementsPerCrossPoint(k,3) = A(j+1,i);
                    obj.elementsPerCrossPoint(k,4) = A(j+1,i+1);
                elseif j == 0
                    obj.crossPointsPerElement(A(1,i),2)   = k;
                    obj.crossPointsPerElement(A(1,i+1),1) = k;

                    obj.elementsPerCrossPoint(k,1) = A(1,i);
                    obj.elementsPerCrossPoint(k,2) = A(1,i+1);
                elseif i == n
                    obj.crossPointsPerElement(A(j,n),4)   = k;
                    obj.crossPointsPerElement(A(j+1,n),2) = k;

                    obj.elementsPerCrossPoint(k,1) = A(j,n);
                    obj.elementsPerCrossPoint(k,2) = A(j+1,n);
                elseif j == n
                    obj.crossPointsPerElement(A(n,i),4)   = k;
                    obj.crossPointsPerElement(A(n,i+1),3) = k;

                    obj.elementsPerCrossPoint(k,1) = A(n,i);
                    obj.elementsPerCrossPoint(k,2) = A(n,i+1);
                elseif i == 0
                    obj.crossPointsPerElement(A(j,1),3)   = k;
                    obj.crossPointsPerElement(A(j+1,1),1) = k;

                    obj.elementsPerCrossPoint(k,1) = A(j,1);
                    obj.elementsPerCrossPoint(k,2) = A(j+1,1);
                end
            end

            obj.enclosedCrossPointsPerMacroFace = {};
            obj.macroFacesPerLevel = zeros(1,0);
        end

        function divide(obj, flag)

            prevGrid = obj.map{end,1};
            A0 = prevGrid{1,1};

            if flag == 0
                if size(A0,1) < 2
                    return
                end
                nextGrid = cell(2*size(prevGrid,1), size(prevGrid,2));
                pointsPerMergedFace = zeros(numel(prevGrid), size(A0,2)+1);
            else
                if size(A0,2) < 2
                    return
                end
                nextGrid = cell(size(prevGrid,1), 2*size(prevGrid,2));
                pointsPerMergedFace = zeros(numel(prevGrid), size(A0,1)+1);
            end

            facePairs = zeros(0,2);
            mergeFace = zeros(0,2);

            if flag == 0
                for c = 1:size(prevGrid,2)
                    for r = 1:size(prevGrid,1)
                        A = prevGrid{r,c};
                        mid = size(A,1)/2;
                        idx = (c-1)*size(prevGrid,1) + r;

                        nextGrid{2*r-1, c} = A(1:mid, :);
                        nextGrid{2*r,   c} = A(mid+1:end, :);

                        for j = 1:size(A,2)
                            facePairs(end+1,:) = [A(mid,j), A(mid+1,j)];
                            mergeFace(end+1,:) = [4,3];
                        end

                        if c == 1
                            pointsPerMergedFace(idx,1) = obj.crossPointsPerElement(A(mid,1),3);
                        end

                        for j = 1:size(A,2)-1
                            pointsPerMergedFace(idx,j+1) = obj.crossPointsPerElement(A(mid,j),4);
                        end

                        if c == size(prevGrid,2)
                            pointsPerMergedFace(idx,size(A,2)+1) = obj.crossPointsPerElement(A(mid,size(A,2)),4);
                        end
                    end
                end
            else
                for r = 1:size(prevGrid,1)
                    for c = 1:size(prevGrid,2)
                        A = prevGrid{r,c};
                        mid = size(A,2)/2;
                        idx = (r-1)*size(prevGrid,2) + c;

                        nextGrid{r, 2*c-1} = A(:, 1:mid);
                        nextGrid{r, 2*c}   = A(:, mid+1:end);

                        for i = 1:size(A,1)
                            facePairs(end+1,:) = [A(i,mid), A(i,mid+1)];
                            mergeFace(end+1,:) = [2,1];
                        end

                        if r == 1
                            pointsPerMergedFace(idx,1) = obj.crossPointsPerElement(A(1,mid),2);
                        end

                        for i = 1:size(A,1)-1
                            pointsPerMergedFace(idx,i+1) = obj.crossPointsPerElement(A(i,mid),4);
                        end

                        if r == size(prevGrid,1)
                            pointsPerMergedFace(idx,size(A,1)+1) = obj.crossPointsPerElement(A(size(A,1),mid),4);
                        end
                    end
                end
            end

            obj.map{end+1,1} = nextGrid;

            obj.divide(1-flag);

            for i = 1:size(facePairs,1)
                if flag == 1
                    obj.facePerElement(facePairs(i,1),2) = size(obj.elementsPerFace,1) + i;
                    obj.facePerElement(facePairs(i,2),1) = size(obj.elementsPerFace,1) + i;
                else
                    obj.facePerElement(facePairs(i,1),4) = size(obj.elementsPerFace,1) + i;
                    obj.facePerElement(facePairs(i,2),3) = size(obj.elementsPerFace,1) + i;
                end
            end

            obj.levels(end+1) = size(facePairs,1);
            obj.elementsPerFace = [obj.elementsPerFace; facePairs];
            obj.elementSidePerFace = [obj.elementSidePerFace; mergeFace];
            obj.enclosedCrossPointsPerMacroFace{1,end+1} = pointsPerMergedFace;
            obj.macroFacesPerLevel(end+1) = size(pointsPerMergedFace,1);
        end

        function calculateReordering(obj, dofsPerFace)

            nFaces = size(obj.elementsPerFace, 1);
            nRawFaceDofs = 2 * nFaces * dofsPerFace;

            pointDofPerElement = zeros(size(obj.crossPointsPerElement));
            pointDofCursor = 0;
            for e = 1:size(obj.crossPointsPerElement,1)
                for c = 1:4
                    if obj.crossPointsPerElement(e,c) ~= 0
                        pointDofCursor = pointDofCursor + 1;
                        pointDofPerElement(e,c) = nRawFaceDofs + pointDofCursor;
                    end
                end
            end

            obj.permutation = zeros(0,1);
            obj.nDofsPerMacroFace = zeros(0,1);

            usedPoint = false(size(pointDofPerElement));

            % 1 = left, 2 = right, 3 = bottom, 4 = top
            % boundary sides count as already numbered
            numberedSide = (obj.facePerElement == 0);

            rawCursor = 0;
            faceId = 0;

            for level = 1:numel(obj.levels)
                facesPerMacroFace = obj.levels(level) / obj.macroFacesPerLevel(level);

                for macroFace = 1:obj.macroFacesPerLevel(level)
                    blockStart = numel(obj.permutation) + 1;

                    for f = 1:facesPerMacroFace
                        faceId = faceId + 1;

                        for q = 1:2
                            e = obj.elementsPerFace(faceId,q);
                            s = obj.elementSidePerFace(faceId,q);

                            % corners:
                            % 1 = bottom-left
                            % 2 = bottom-right
                            % 3 = top-left
                            % 4 = top-right

                            if s == 1
                                cStart = 1;
                                cEnd = 3;
                            elseif s == 2
                                cStart = 2;
                                cEnd = 4;
                            elseif s == 3
                                cStart = 1;
                                cEnd = 2;
                            else
                                cStart = 3;
                                cEnd = 4;
                            end

                            if cStart == 1
                                otherStartSide = 3;
                            elseif cStart == 2
                                otherStartSide = 3;
                            elseif cStart == 3
                                otherStartSide = 4;
                            else
                                otherStartSide = 4;
                            end
                            if otherStartSide == s
                                if cStart == 1 || cStart == 3
                                    otherStartSide = 1;
                                else
                                    otherStartSide = 2;
                                end
                            end

                            if cEnd == 1
                                otherEndSide = 3;
                            elseif cEnd == 2
                                otherEndSide = 3;
                            elseif cEnd == 3
                                otherEndSide = 4;
                            else
                                otherEndSide = 4;
                            end
                            if otherEndSide == s
                                if cEnd == 1 || cEnd == 3
                                    otherEndSide = 1;
                                else
                                    otherEndSide = 2;
                                end
                            end

                            id = pointDofPerElement(e,cStart);
                            if id ~= 0 && ~usedPoint(e,cStart) && numberedSide(e,otherStartSide)
                                obj.permutation(end+1,1) = id;
                                usedPoint(e,cStart) = true;
                            end

                            ids = rawCursor + (1:dofsPerFace);
                            obj.permutation = [obj.permutation; ids(:)];
                            rawCursor = rawCursor + dofsPerFace;

                            numberedSide(e,s) = true;

                            id = pointDofPerElement(e,cEnd);
                            if id ~= 0 && ~usedPoint(e,cEnd) && numberedSide(e,otherEndSide)
                                obj.permutation(end+1,1) = id;
                                usedPoint(e,cEnd) = true;
                            end
                        end
                    end

                    obj.nDofsPerMacroFace(end+1,1) = numel(obj.permutation) - blockStart + 1;
                end
            end

            for e = 1:size(pointDofPerElement,1)
                for c = 1:4
                    id = pointDofPerElement(e,c);
                    if id ~= 0 && ~usedPoint(e,c)
                        obj.permutation(end+1,1) = id;
                        usedPoint(e,c) = true;
                    end
                end
            end

        end

    end
end