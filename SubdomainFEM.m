classdef SubdomainFEM
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
        eta
        poincareSteklovOperator
    end

    methods
        function obj = SubdomainFEM(div, f_fun, ax, bx, ay, by, eta, c0, poincareSteklovOperator)
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

            if ax >= bx
                error('Subdomain: ax must be < bx.');
            end
            if ay >= by
                error('Subdomain: ay must be < by.');
            end

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

            if obj.poincareSteklovOperator ~= "DtN"
                error('Subdomain: ItI requested, but this FEM Schur constructor currently implements the DtN-style Schur complement only.');
            end

            [px,py,K,rhs0] = get_fem(div, c0fun, f_fun, ax, bx, ay, by, 2^div, 2^div);

            obj.px = px;
            obj.py = py;
            obj.K = full(K);
            obj.rhs0 = rhs0;

            tol = 1.0e-14;
            obj.idx_interior = find(~((abs(px-ax)<tol) | (abs(px-bx)<tol) | (abs(py-ay)<tol) | (abs(py-by)<tol)));
            obj.idx_corners  = find(((abs(px-ax)<tol) | (abs(px-bx)<tol)) & ((abs(py-ay)<tol) | (abs(py-by)<tol)));
            obj.idx_left     = find((abs(px-ax)<tol) & ~((abs(py-ay)<tol) | (abs(py-by)<tol)));
            obj.idx_right    = find((abs(px-bx)<tol) & ~((abs(py-ay)<tol) | (abs(py-by)<tol)));
            obj.idx_bottom   = find((abs(py-ay)<tol) & ~((abs(px-ax)<tol) | (abs(px-bx)<tol)));
            obj.idx_top      = find((abs(py-by)<tol) & ~((abs(px-ax)<tol) | (abs(px-bx)<tol)));
            obj.idx_boundary = [obj.idx_left; obj.idx_right; obj.idx_bottom; obj.idx_top; obj.idx_corners];

            obj.b = [numel(obj.idx_left), numel(obj.idx_right), numel(obj.idx_bottom), numel(obj.idx_top), numel(obj.idx_corners)];

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