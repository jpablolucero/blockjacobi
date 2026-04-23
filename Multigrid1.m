classdef Multigrid1 < handle
    properties
        A
        invA
        B
        C
        D
        M
        separatorBlockSizes
        separatorPermutation
        tol
        maxit
        restart
        m
    end

    methods
        function obj = Multigrid1(A_in, nd, lmax, tol, maxit, restart, m)
            obj.separatorBlockSizes = nd.separatorBlockSizes;
            obj.separatorPermutation = nd.separatorPermutation;
            obj.tol = tol;
            obj.maxit = maxit;
            obj.restart = restart;
            obj.m = m;

            obj.A = {};
            obj.invA = {};
            obj.B = {};
            obj.C = {};
            obj.D = {};
            obj.M = {};

            p = obj.separatorPermutation{1};
            obj.M{1} = sparse(A_in(p, p));

            ilevel = 1;

            while ilevel < lmax
                nSeparator = sum(obj.separatorBlockSizes{ilevel});

                A = sparse(obj.M{ilevel}(1:nSeparator, 1:nSeparator));
                B = sparse(obj.M{ilevel}(1:nSeparator, nSeparator + 1:end));
                C = sparse(obj.M{ilevel}(nSeparator + 1:end, 1:nSeparator));
                D = sparse(obj.M{ilevel}(nSeparator + 1:end, nSeparator + 1:end));

                S = sparse(D - C * (A \ B));

                obj.A{ilevel} = A;
                obj.invA{ilevel} = inv(A);
                obj.B{ilevel} = B;
                obj.C{ilevel} = C;
                obj.D{ilevel} = D;

                if ilevel + 1 <= numel(obj.separatorPermutation)
                    p = obj.separatorPermutation{ilevel + 1};
                    obj.M{ilevel + 1} = sparse(S(p, p));
                else
                    obj.M{ilevel + 1} = sparse(S);
                end

                ilevel = ilevel + 1;
            end
        end

        function z = P(obj, y, level)
            if level + 1 <= numel(obj.separatorPermutation)
                p = obj.separatorPermutation{level + 1};
                yUnpermuted = zeros(size(y));
                yUnpermuted(p) = y;
            else
                yUnpermuted = y;
            end

            z = [-(obj.invA{level} * (obj.B{level} * yUnpermuted)); yUnpermuted];
        end

        function y = R(obj, x, level)
            nSeparator = size(obj.B{level}, 1);

            x1 = x(1:nSeparator);
            x2 = x(nSeparator + 1:end);

            yUnpermuted = -obj.C{level} * (obj.invA{level} * x1) + x2;

            if level + 1 <= numel(obj.separatorPermutation)
                p = obj.separatorPermutation{level + 1};
                y = yUnpermuted(p);
            else
                y = yUnpermuted;
            end
        end

        function x = vcycle(obj, g, level)
            if nargin < 3
                p = obj.separatorPermutation{1};
                x0 = obj.vcycle(g(p), 1);
                x = zeros(size(x0));
                x(p) = x0;
                return
            end

            if level < length(obj.M)
                nSeparator = size(obj.B{level}, 1);

                x = zeros(size(g));

                r = g - obj.M{level} * x;
                x = x + [obj.invA{level} * r(1:nSeparator); zeros(length(g) - nSeparator, 1)];

                r = obj.R(g - obj.M{level} * x, level);
                e = obj.vcycle(r, level + 1);
                x = x + obj.P(e, level);
            else
                x = obj.M{level} \ g;
            end
        end
    end
end