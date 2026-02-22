classdef Subdomain
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
        function obj = Subdomain(div, f_fun, ax, bx, ay, by, tol, mx, my)
            if nargin < 10, my = 2^div; end
            if nargin < 9,  mx = 2^div; end
            if nargin < 8,  tol = 1e-14; end

            [~,~,~,px,py,K,rhs0] = get_fem(div, @(x,y) 0*x, f_fun, ax, bx, ay, by, tol, mx, my);
            obj.px = px;
            obj.py = py;
            obj.K   = full(K);
            obj.rhs0 = rhs0;

            obj.idx_interior = find(~((abs(px-ax)<tol)|(abs(px-bx)<tol)|(abs(py-ay)<tol)|(abs(py-by)<tol)));
            obj.idx_corners  = find(((abs(px-ax)<tol)|(abs(px-bx)<tol)) & ((abs(py-ay)<tol)|(abs(py-by)<tol)));
            obj.idx_left     = find((abs(px-ax)<tol) & ~((abs(py-ay)<tol)|(abs(py-by)<tol)));
            obj.idx_right    = find((abs(px-bx)<tol) & ~((abs(py-ay)<tol)|(abs(py-by)<tol)));
            obj.idx_bottom   = find((abs(py-ay)<tol) & ~((abs(px-ax)<tol)|(abs(px-bx)<tol)));
            obj.idx_top      = find((abs(py-by)<tol) & ~((abs(px-ax)<tol)|(abs(px-bx)<tol)));
            obj.idx_boundary = [obj.idx_left; obj.idx_right; obj.idx_bottom; obj.idx_top; obj.idx_corners];

            obj.b = [numel(obj.idx_left), numel(obj.idx_right), numel(obj.idx_bottom), numel(obj.idx_top), numel(obj.idx_corners)];

            obj.A = obj.K(obj.idx_interior, obj.idx_interior);
            obj.B = obj.K(obj.idx_interior, obj.idx_boundary);
            obj.C = obj.K(obj.idx_boundary, obj.idx_interior);
            obj.D = obj.K(obj.idx_boundary, obj.idx_boundary);

            obj.S   = obj.D - obj.C * (obj.A \ obj.B);
            obj.T = mat2cell(obj.S, obj.b, obj.b);

            obj.rhs_interior = rhs0(obj.idx_interior);
            h_full = rhs0(obj.idx_boundary) - obj.C * (obj.A \ obj.rhs_interior);
            obj.h = mat2cell(h_full, obj.b, 1);
        end
        function indices = idx(obj,block)
            if block==1
                indices = obj.idx_left;
            elseif block == 2
                indices = obj.idx_right;
            elseif block == 3
                indices = obj.idx_bottom;
            elseif block == 4
                indices = obj.idx_top;
            elseif block == 5
                indices = obj.idx_corners;
            end
        end
    end
end
