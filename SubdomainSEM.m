classdef SubdomainSEM
    properties (SetAccess=public)
        idx_interior
        idx_left
        idx_right
        idx_bottom
        idx_top
        idx_corners
        idx_boundary
        b
        A
        B
        C
        D
        S
        T
        K
        rhs0
        rhs_interior
        h
        px
        py
    end

    methods
        function obj = SubdomainSEM(div, f_fun, ax, bx, ay, by, mx, my)
            if nargin < 6
                error('SubdomainSEM: You must specify div, f_fun, ax, bx, ay, by.');
            end
            if nargin < 7
                mx = 1;
            end
            if nargin < 8
                my = 1;
            end

            if ax >= bx
                error('SubdomainSEM: ax must be < bx.');
            end
            if ay >= by
                error('SubdomainSEM: ay must be < by.');
            end
            if mx ~= 1
                error('SubdomainSEM expects one element per subdomain: mx=1.');
            end
            if my ~= 1
                error('SubdomainSEM expects one element per subdomain: my=1.');
            end

            p = 2^div;
            quad = Quadrature(p);

            x = quad.Points;
            y = quad.Points;

            xscale = (bx - ax) / 2;
            yscale = (by - ay) / 2;

            xCoords = xscale * x + (bx + ax) / 2;
            yCoords = yscale * y + (by + ay) / 2;

            [X, Y] = meshgrid(xCoords, yCoords);
            obj.px = reshape(X.', [], 1);
            obj.py = reshape(Y.', [], 1);

            wx = quad.Weights;
            wy = quad.Weights;

            xWeights = xscale * wx;
            yWeights = yscale * wy;

            Dx = SubdomainSEM.computeDifferentiationMatrix(x, quad.Polynomials) / xscale;
            Dy = SubdomainSEM.computeDifferentiationMatrix(y, quad.Polynomials) / yscale;

            Mx = spdiags(xWeights, 0, numel(xWeights), numel(xWeights));
            My = spdiags(yWeights, 0, numel(yWeights), numel(yWeights));

            obj.K = full(kron(Dx' * Mx * Dx, My) + kron(Mx, Dy' * My * Dy));

            rhs_vals = f_fun(obj.px, obj.py);
            obj.rhs0 = kron(xWeights, yWeights) .* rhs_vals;

            tol = 1.0e-14;

            on_left = abs(obj.px - ax) < tol;
            on_right = abs(obj.px - bx) < tol;
            on_bottom = abs(obj.py - ay) < tol;
            on_top = abs(obj.py - by) < tol;

            obj.idx_interior = find(~(on_left | on_right | on_bottom | on_top));
            obj.idx_corners = find((on_left | on_right) & (on_bottom | on_top));

            obj.idx_left = find(on_left & ~(on_bottom | on_top));
            obj.idx_right = find(on_right & ~(on_bottom | on_top));
            obj.idx_bottom = find(on_bottom & ~(on_left | on_right));
            obj.idx_top = find(on_top & ~(on_left | on_right));

            obj.idx_boundary = [obj.idx_left;
                                obj.idx_right;
                                obj.idx_bottom;
                                obj.idx_top;
                                obj.idx_corners];

            obj.b = [numel(obj.idx_left);
                     numel(obj.idx_right);
                     numel(obj.idx_bottom);
                     numel(obj.idx_top);
                     numel(obj.idx_corners)];

            obj.A = obj.K(obj.idx_interior, obj.idx_interior);
            obj.B = obj.K(obj.idx_interior, obj.idx_boundary);
            obj.C = obj.K(obj.idx_boundary, obj.idx_interior);
            obj.D = obj.K(obj.idx_boundary, obj.idx_boundary);

            obj.S = obj.D - obj.C * (obj.A \ obj.B);
            obj.T = mat2cell(obj.S, obj.b, obj.b);

            obj.rhs_interior = obj.rhs0(obj.idx_interior);

            h_full = obj.rhs0(obj.idx_boundary) - obj.C * (obj.A \ obj.rhs_interior);
            obj.h = mat2cell(h_full, obj.b, 1);
        end

        function indices = idx(obj, block)
            if block == 1
                indices = obj.idx_left;
                return;
            end
            if block == 2
                indices = obj.idx_right;
                return;
            end
            if block == 3
                indices = obj.idx_bottom;
                return;
            end
            if block == 4
                indices = obj.idx_top;
                return;
            end
            if block == 5
                indices = obj.idx_corners;
                return;
            end
            error('SubdomainSEM.idx: block must be 1..5.');
        end
    end

    methods (Static)
        function D = computeDifferentiationMatrix(x, P)
            N1 = length(x);

            X = repmat(x, 1, N1);
            Xdiff = X - X' + eye(N1);

            L = repmat(P(:, end), 1, N1);
            L(1:(N1+1):N1*N1) = 1;

            D = (L ./ (Xdiff .* L'));
            D(1:(N1+1):N1*N1) = 0;

            N = N1 - 1;
            D(1, 1) = -(N1 * N) / 4;
            D(end, end) = (N1 * N) / 4;
        end
    end
end