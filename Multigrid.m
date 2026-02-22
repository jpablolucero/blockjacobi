classdef Multigrid < handle
    properties
        N
        A
        B
        C
        D
        M
        bSize
        tol
        m
    end
    methods
        function obj = Multigrid(A_in,nd,lmax,tol,m)
            obj.bSize = nd.nDofsPerMacroFace(:);
            obj.m = m;
            obj.A     = {};
            obj.M     = {};
            obj.B     = {};
            obj.C     = {};
            obj.D     = {};
            M         = A_in;
            obj.M{1}  = M;
            ilevel    = 1;
            while ilevel < lmax
                n_b    = nd.macroFacesPerLevel(ilevel);
                bL     = obj.bSize(ilevel);
                ilevel = ilevel + 1;

                A = obj.M{end}(1:n_b*bL,           1:n_b*bL);
                B = obj.M{end}(1:n_b*bL,           n_b*bL+1:end);
                C = obj.M{end}(n_b*bL+1:end,       1:n_b*bL);
                D = obj.M{end}(n_b*bL+1:end,       n_b*bL+1:end);

                obj.M{end+1} = sparse(D - C*(A\B));
                obj.A{end+1} = A;
                obj.B{end+1} = B;
                obj.C{end+1} = C;
                obj.D{end+1} = D;
            end
        end

        function z = P(obj, y, level)
            z  = [-(obj.A{level} \ (obj.B{level} * y)); y];
        end

        function y = R(obj, x, level)
            nA = size(obj.B{level},1);
            x1 = x(1:nA);
            x2 = x(nA+1:end);
            y  = -obj.C{level} * (obj.A{level} \ x1) + x2;
        end

        function x = vcycle(obj,g,level)
            if nargin < 3
                level = 1;
            end

            rlx = obj.OptimalRelaxations(obj.m);

            if (level<length(obj.M))
                
                x = zeros(size(g));

                for i = 1:length(rlx)
                    nA = size(obj.B{level},1);
                    r  = g - obj.M{level}*x;
                    x  = x + rlx(i) * [ obj.M{level}(1:nA,1:nA)\r(1:nA); obj.M{level}(nA+1:end,nA+1:end)\r(nA+1:end) ];
                end

                r    = obj.R(g - obj.M{level}*x, level);

                Minv = @(z) obj.vcycle(z, level+1);
                rho  = 1 - 1/(2*obj.m + 1)^2;

                MinvR   = Minv(r);
                e    = (1 + 1/rho) * MinvR - (1/rho) * Minv(obj.M{level+1} * MinvR);

                x = x + obj.P(e, level);

                for i = 1:length(rlx)
                    nA = size(obj.B{level},1);
                    r  = g - obj.M{level}*x;
                    x  = x + rlx(length(rlx)+1-i) * [ obj.M{level}(1:nA,1:nA)\r(1:nA); obj.M{level}(nA+1:end,nA+1:end)\r(nA+1:end) ];
                end
           
            else       
                x = obj.M{level}\g;
            end
        end

        function rlx = OptimalRelaxations(~,m)
            j = 1:m;
            rlx = 1 ./ (1 - cos(2*pi*j./(2*m + 1)));
        end
    end
end