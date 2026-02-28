classdef Subdomain < handle
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
        ax
        ay
        bx
        by
    end

    methods
        function obj = Subdomain(div, rhs_fun, ax, bx, ay, by, eta, c0, poincareSteklovOperator, method)
            if nargin < 6
                error('Subdomain: You must specify div, f_fun, ax, bx, ay, by.');
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
                error('Subdomain: ax must be < bx.');
            end
            if ay >= by
                error('Subdomain: ay must be < by.');
            end

            obj.ax = ax;
            obj.ay = ay;
            obj.bx = bx;
            obj.by = by;

            if (poincareSteklovOperator == "DtN")
                obj.poincareSteklovOperator = "DtN";
            elseif (poincareSteklovOperator == "ItI")
                obj.poincareSteklovOperator = "ItI";
            else
                error('Subdomain: Unknown Poincare Steklov Operator.');
            end
            obj.eta = eta;

            if isa(c0, 'function_handle')
                c0fun = c0;
            else
                c0const = c0;
                c0fun = @(x,y) c0const + 0 .* x;
            end

            [obj.px,obj.py,obj.K,obj.rhs,~] = method(div, c0fun, rhs_fun, ax, bx, ay, by);

            tol = 1.0e-14;
            obj.idx_interior = find(~((abs(obj.px-ax)<tol) |   (abs(obj.px-bx)<tol)  |  (abs(obj.py-ay)<tol) | (abs(obj.py-by)<tol)));
            obj.idx_corners  = find( ((abs(obj.px-ax)<tol) |   (abs(obj.px-bx)<tol)) & ((abs(obj.py-ay)<tol) | (abs(obj.py-by)<tol)));
            obj.idx_left     = find(  (abs(obj.px-ax)<tol) & ~((abs(obj.py-ay)<tol)  |  (abs(obj.py-by)<tol)));
            obj.idx_right    = find(  (abs(obj.px-bx)<tol) & ~((abs(obj.py-ay)<tol)  |  (abs(obj.py-by)<tol)));
            obj.idx_bottom   = find(  (abs(obj.py-ay)<tol) & ~((abs(obj.px-ax)<tol)  |  (abs(obj.px-bx)<tol)));
            obj.idx_top      = find(  (abs(obj.py-by)<tol) & ~((abs(obj.px-ax)<tol)  |  (abs(obj.px-bx)<tol)));
            obj.idx_boundary = [obj.idx_left; obj.idx_right; obj.idx_bottom; obj.idx_top; obj.idx_corners];

            obj.b = [numel(obj.idx_left), numel(obj.idx_right), numel(obj.idx_bottom), numel(obj.idx_top), numel(obj.idx_corners)];

            obj.A = obj.K(obj.idx_interior, obj.idx_interior);
            obj.B = obj.K(obj.idx_interior, obj.idx_boundary);
            obj.C = obj.K(obj.idx_boundary, obj.idx_interior);
            obj.D = obj.K(obj.idx_boundary, obj.idx_boundary);

            obj.S = obj.D - obj.C * (obj.A \ obj.B);

            h_full = obj.rhs(obj.idx_boundary) - obj.C * (obj.A \ obj.rhs(obj.idx_interior));

            if obj.poincareSteklovOperator == "DtN"
                obj.T = mat2cell(obj.S, obj.b, obj.b);
                obj.h = mat2cell(h_full, obj.b, 1);
            else
                error('Subdomain: Unknown Poincare Steklov Operator.');
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

                    val = 0;

                    if (abs(obj.px(id) - 0) < tol) && (abs(obj.py(id) - 0) < tol)
                        val = BC{3}(obj.px(id), obj.py(id));
                    elseif (abs(obj.px(id) - 1) < tol) && (abs(obj.py(id) - 0) < tol)
                        val = BC{2}(obj.px(id), obj.py(id));
                    elseif (abs(obj.px(id) - 1) < tol) && (abs(obj.py(id) - 1) < tol)
                        val = BC{4}(obj.px(id), obj.py(id));
                    elseif (abs(obj.px(id) - 0) < tol) && (abs(obj.py(id) - 1) < tol)
                        val = BC{1}(obj.px(id), obj.py(id));
                    elseif abs(obj.px(id) - 0) < tol
                        val = BC{1}(obj.px(id), obj.py(id));
                    elseif abs(obj.px(id) - 1) < tol
                        val = BC{2}(obj.px(id), obj.py(id));
                    elseif abs(obj.py(id) - 0) < tol
                        val = BC{3}(obj.px(id), obj.py(id));
                    elseif abs(obj.py(id) - 1) < tol
                        val = BC{4}(obj.px(id), obj.py(id));
                    end

                    for r = 1:5
                        obj.h{r} = obj.h{r} - obj.T{r,5}(:,c) * val;
                    end
                end
            end
        end

        function indices = idx(obj,block)
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
            error('Subdomain.idx: block must be 1..5.');
        end
    end
end