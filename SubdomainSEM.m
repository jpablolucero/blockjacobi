classdef SubdomainSEM < handle
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
        rhs
        h
        px
        py
        eta
        poincareSteklovOperator
    end

    methods
        function obj = SubdomainSEM(div, f_fun, ax, bx, ay, by, eta, c0, poincareSteklovOperator, BC)
            if nargin < 6
                error('SubdomainSEM: You must specify div, f_fun, ax, bx, ay, by.');
            end
            if nargin < 7
                eta = 0;
            end
            if nargin < 8
                c0 = 0;
            end
            if nargin < 9
                poincareSteklovOperator = "DtN";
            end
            if nargin < 10
                BC = [];
            end

            if ax >= bx
                error('SubdomainSEM: ax must be < bx.');
            end
            if ay >= by
                error('SubdomainSEM: ay must be < by.');
            end

            if (poincareSteklovOperator == "DtN")
                obj.poincareSteklovOperator = "DtN";
            elseif (poincareSteklovOperator == "ItI")
                obj.poincareSteklovOperator = "ItI";
            else
                error('SubdomainSEM: Unknown Poincare Steklov Operator.');
            end
            obj.eta = eta;

            if isa(c0, 'function_handle')
                c0fun = c0;
            else
                c0const = c0;
                c0fun = @(x,y) c0const + 0 .* x;
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

            c0vals = c0fun(obj.px, obj.py);

            obj.K = full( ...
                kron(Dx' * Mx * Dx, My) + ...
                kron(Mx, Dy' * My * Dy) + ...
                c0vals .* kron(My, Mx) );

            rhs_vals = f_fun(obj.px, obj.py);
            obj.rhs = kron(xWeights, yWeights) .* rhs_vals;

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

            if obj.poincareSteklovOperator ~= "DtN"
                error('SubdomainSEM: ItI requested, but this SEM Schur constructor currently implements the DtN-style Schur complement only.');
            end

            obj.A = obj.K(obj.idx_interior, obj.idx_interior);
            obj.B = obj.K(obj.idx_interior, obj.idx_boundary);
            obj.C = obj.K(obj.idx_boundary, obj.idx_interior);
            obj.D = obj.K(obj.idx_boundary, obj.idx_boundary);

            obj.S = obj.D - obj.C * (obj.A \ obj.B);
            obj.T = mat2cell(obj.S, obj.b, obj.b);

            h_full = obj.rhs(obj.idx_boundary) - obj.C * (obj.A \ obj.rhs(obj.idx_interior));
            obj.h = mat2cell(h_full, obj.b, 1);

            if ~isempty(BC)
                obj.setBoundaryCondition(BC, tol);
            end
        end

        function setBoundaryCondition(obj, BC, tol)
            if nargin < 3
                tol = 1.0e-12;
            end

            if any(abs(obj.px(obj.idx_left) - 0) < tol)
                for r = 1:5
                    obj.h{r} = obj.h{r} - obj.T{r,1} * BC{1}(obj.px(obj.idx_left), obj.py(obj.idx_left));
                end
            end

            if any(abs(obj.px(obj.idx_right) - 1) < tol)
                for r = 1:5
                    obj.h{r} = obj.h{r} - obj.T{r,2} * BC{2}(obj.px(obj.idx_right), obj.py(obj.idx_right));
                end
            end

            if any(abs(obj.py(obj.idx_bottom) - 0) < tol)
                for r = 1:5
                    obj.h{r} = obj.h{r} - obj.T{r,3} * BC{3}(obj.px(obj.idx_bottom), obj.py(obj.idx_bottom));
                end
            end

            if any(abs(obj.py(obj.idx_top) - 1) < tol)
                for r = 1:5
                    obj.h{r} = obj.h{r} - obj.T{r,4} * BC{4}(obj.px(obj.idx_top), obj.py(obj.idx_top));
                end
            end

            for c = 1:4
                id = obj.idx_corners(c);

                if abs(obj.px(id) - 0) < tol || abs(obj.px(id) - 1) < tol || abs(obj.py(id) - 0) < tol || abs(obj.py(id) - 1) < tol
                    vals = [];

                    if abs(obj.px(id) - 0) < tol
                        vals(end+1,1) = BC{1}(obj.px(id), obj.py(id));
                    end

                    if abs(obj.px(id) - 1) < tol
                        vals(end+1,1) = BC{2}(obj.px(id), obj.py(id));
                    end

                    if abs(obj.py(id) - 0) < tol
                        vals(end+1,1) = BC{3}(obj.px(id), obj.py(id));
                    end

                    if abs(obj.py(id) - 1) < tol
                        vals(end+1,1) = BC{4}(obj.px(id), obj.py(id));
                    end

                    for r = 1:5
                        obj.h{r} = obj.h{r} - obj.T{r,5}(:,c) * mean(vals);
                    end
                end
            end
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