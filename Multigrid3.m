classdef Multigrid3 < handle
    properties
        N
        A
        B
        C
        D
        M
        bSize
        tol
        maxit
        restart
        m
    end
    methods
        function obj = Multigrid3(A_in,nd,lmax,tol,maxit,restart,m)
            obj.bSize = nd.nDofsPerMacroFace;
            obj.tol = tol; obj.maxit = maxit ; obj.restart = restart;obj.m = m;
            obj.A = {};
            obj.M = {};
            obj.B = {};
            obj.C = {};
            obj.D = {};

            obj.M{1} = sparse(A_in);
            ilevel = 1;

            while ilevel < lmax
                nA = sum(obj.bSize{ilevel});

                A = sparse(obj.M{end}(1:nA,     1:nA));
                B = sparse(obj.M{end}(1:nA,     nA+1:end));
                C = sparse(obj.M{end}(nA+1:end, 1:nA));
                D = sparse(obj.M{end}(nA+1:end, nA+1:end));

                obj.M{end+1} = sparse(D - C * (A \ B));
                obj.A{end+1} = A;
                obj.B{end+1} = B;
                obj.C{end+1} = C;
                obj.D{end+1} = D;

                ilevel = ilevel + 1;
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
            rho  = 1 - 1/(2*obj.m + 1)^2;

            if (level<length(obj.M))
                
                x = zeros(size(g));

                for i = 1:length(rlx)
                    nA = size(obj.B{level},1);
                    r  = g - obj.M{level}*x;
                    x  = x + rlx(i) * [ obj.M{level}(1:nA,1:nA)\r(1:nA); obj.M{level}(nA+1:end,nA+1:end)\r(nA+1:end) ];
                end

                r    = obj.R(g - obj.M{level}*x, level);

                MinvR = obj.vcycle(r, level+1);
                e    = (1 + 1/rho) * MinvR - (1/rho) * obj.vcycle(obj.M{level+1} * MinvR, level+1);
               
                x = x + obj.P(e, level);

                for i = 1:length(rlx)
                    nA = size(obj.B{level},1);
                    r  = g - obj.M{level}*x;
                    x  = x + rlx(length(rlx)+1-i) * [ obj.M{level}(1:nA,1:nA)\r(1:nA); obj.M{level}(nA+1:end,nA+1:end)\r(nA+1:end) ];
                end
           
            else       
                % x = obj.M{level} \ g;

                [x, ~, ~] = gmres( ...
                    obj.M{level}, ...
                    g, ...
                    obj.restart, ...
                    obj.tol, ...
                    obj.maxit, ...
                    [], ...
                    [], ...
                    zeros(size(g)));

            end
        end

        function rlx = OptimalRelaxations(~,m)
            j = 1:m;
            rlx = 1 ./ (1 - cos(2*pi*j./(2*m + 1)));
        end
    end
end