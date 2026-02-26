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
    end

    methods
        function obj = NestedDissection(arg)
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

            obj.map = cell(1,1);
            obj.map{1,1} = cell(1,1);
            obj.map{1,1}{1,1} = A;

            obj.elementsPerFace = zeros(0,2);
            obj.elementSidePerFace = zeros(0,2);

            obj.facePerElement = zeros(2^obj.N*2^obj.N,4);

            obj.levels = [];

            obj.crossPointGrid = [];
            for j=1:2^obj.N-1
                for i=1:2^obj.N-1
                    obj.crossPointGrid(end+1,:) = [i,j];
                end
            end

            obj.crossPointsPerElement = zeros(2^obj.N*2^obj.N,4);
            obj.elementsPerCrossPoint = zeros((2^obj.N-1)*(2^obj.N-1),4);

            for i=1:(2^obj.N-1)*(2^obj.N-1)
                obj.crossPointsPerElement(A(obj.crossPointGrid(i,2),obj.crossPointGrid(i,1)),4) = i;
                obj.crossPointsPerElement(A(obj.crossPointGrid(i,2),obj.crossPointGrid(i,1)+1),3) = i;
                obj.crossPointsPerElement(A(obj.crossPointGrid(i,2)+1,obj.crossPointGrid(i,1)),2) = i;
                obj.crossPointsPerElement(A(obj.crossPointGrid(i,2)+1,obj.crossPointGrid(i,1)+1),1) = i;

                obj.elementsPerCrossPoint(i,1) = A(obj.crossPointGrid(i,2),obj.crossPointGrid(i,1));
                obj.elementsPerCrossPoint(i,2) = A(obj.crossPointGrid(i,2),obj.crossPointGrid(i,1)+1);
                obj.elementsPerCrossPoint(i,3) = A(obj.crossPointGrid(i,2)+1,obj.crossPointGrid(i,1));
                obj.elementsPerCrossPoint(i,4) = A(obj.crossPointGrid(i,2)+1,obj.crossPointGrid(i,1)+1);
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
            else
                if size(A0,2) < 2
                    return
                end
                nextGrid = cell(size(prevGrid,1), 2*size(prevGrid,2));
            end

            facePairs = zeros(0,2);
            mergeFace = zeros(0,2);

            if flag == 0
                verticesPerMergedFace = zeros(numel(prevGrid), size(A0,2)-1);
            else
                verticesPerMergedFace = zeros(numel(prevGrid), size(A0,1)-1);
            end

            if flag == 0
                for c = 1:size(prevGrid,2)
                    for r = 1:size(prevGrid,1)
                        A = prevGrid{r,c};
                        mid = size(A,1)/2;
                        nextGrid{2*r-1, c} = A(1:mid, :);
                        nextGrid{2*r,   c} = A(mid+1:end, :);
    
                        for j = 1:size(A,2)
                            facePairs(end+1,:) = [A(mid,j), A(mid+1,j)];   % bottom, top (since row is -y)
                            mergeFace(end+1,:) = [4,3];
                        end

                        for j = 1:size(A,2)-1
                            verticesPerMergedFace((c-1)*size(prevGrid,1) + r,j) = obj.crossPointsPerElement(A(mid,j),4);
                        end

                    end
                end
            else
                for r = 1:size(prevGrid,1)
                    for c = 1:size(prevGrid,2)
                        A = prevGrid{r,c};
                        mid = size(A,2)/2;
                        nextGrid{r, 2*c-1} = A(:, 1:mid);
                        nextGrid{r, 2*c}   = A(:, mid+1:end);

                        for i = 1:size(A,1)
                            facePairs(end+1,:) = [A(i,mid), A(i,mid+1)];   % left, right
                            mergeFace(end+1,:) = [2,1];
                        end

                        for i = 1:size(A,1)-1
                            verticesPerMergedFace((r-1)*size(prevGrid,2) + c,i) = obj.crossPointsPerElement(A(i,mid),4);
                        end

                    end
                end
            end

            obj.map{end+1,1} = nextGrid;

            obj.divide(1-flag);

            for i=1:size(facePairs,1)
                if flag==1
                    obj.facePerElement(facePairs(i,1),2) = size(obj.elementsPerFace,1) + i;
                    obj.facePerElement(facePairs(i,2),1) = size(obj.elementsPerFace,1) + i;
                else
                    obj.facePerElement(facePairs(i,1),4) = size(obj.elementsPerFace,1) + i ;
                    obj.facePerElement(facePairs(i,2),3) = size(obj.elementsPerFace,1) + i ;
                end
            end

            obj.levels(end+1) = size(facePairs,1);
            obj.elementsPerFace = [obj.elementsPerFace; facePairs];
            obj.elementSidePerFace = [obj.elementSidePerFace; mergeFace];

            obj.enclosedCrossPointsPerMacroFace{1,end+1} = verticesPerMergedFace;

            obj.macroFacesPerLevel(end+1) = size(verticesPerMergedFace,1);

        end

        function calculateReordering(obj, dofsPerFace)

            nFaceDofsTotal = size(obj.elementsPerFace, 1) * dofsPerFace;
            
            obj.permutation = [];
            obj.nDofsPerMacroFace = [];

            faceDofCursor  = 0;
            for level = 1:numel(obj.levels)
                facesPerMacroFace        = obj.levels(level) / obj.macroFacesPerLevel(level);
                nCrossPointsPerMacroFace = size(obj.enclosedCrossPointsPerMacroFace{level}, 2);
                
                obj.nDofsPerMacroFace(end+1,1)   = facesPerMacroFace * dofsPerFace + nCrossPointsPerMacroFace;

                for macroFace = 1:obj.macroFacesPerLevel(level)

                    nFaces = facesPerMacroFace;
                    cps    = obj.enclosedCrossPointsPerMacroFace{level}(macroFace, :);
                    nCps   = numel(cps);

                    for f = 1:nFaces
                        faceDofs = (faceDofCursor + (f-1)*dofsPerFace + 1):(faceDofCursor + f*dofsPerFace);
                        obj.permutation = [obj.permutation; faceDofs(:)];

                        if f <= nCps
                            crossPointDof = nFaceDofsTotal + cps(f);
                            obj.permutation = [obj.permutation; crossPointDof];
                        end
                    end

                    faceDofCursor = faceDofCursor + nFaces * dofsPerFace;
                end
            
            end

        end

    end
end
