clear;

k   = 2;
eta = k;
div = 5;

[u_exact,rhs,c0,BC,~,IBC,~,~] = DtNTest1(k);

ax = 0; bx = 1;
ay = 0; by = 1;

sd = Subdomain(div, rhs, ax, bx, ay, by, 0, c0, "DtN", @get_sem);

T  = cell2mat(sd.T(1:5,1:5));
h  = cell2mat(sd.h(1:5));
Mb = sd.Mb;

S = Mb \ (T - 1i*eta*Mb) * ((T + 1i*eta*Mb) \ Mb);
c = -(2i*eta) * ((T + 1i*eta*Mb) \ h);

ii = cell(5,1);

ii{1} = IBC{1}(sd.px(sd.idx_left),   sd.py(sd.idx_left));
ii{2} = IBC{2}(sd.px(sd.idx_right),  sd.py(sd.idx_right));
ii{3} = IBC{3}(sd.px(sd.idx_bottom), sd.py(sd.idx_bottom));
ii{4} = IBC{4}(sd.px(sd.idx_top),    sd.py(sd.idx_top));

x1 = sd.px(sd.idx_corners(1)); y1 = sd.py(sd.idx_corners(1));
x2 = sd.px(sd.idx_corners(2)); y2 = sd.py(sd.idx_corners(2));
x3 = sd.px(sd.idx_corners(3)); y3 = sd.py(sd.idx_corners(3));
x4 = sd.px(sd.idx_corners(4)); y4 = sd.py(sd.idx_corners(4));

ii{5} = [0.5*(IBC{1}(x1,y1) + IBC{3}(x1,y1));
         0.5*(IBC{2}(x2,y2) + IBC{3}(x2,y2));
         0.5*(IBC{1}(x3,y3) + IBC{4}(x3,y3));
         0.5*(IBC{2}(x4,y4) + IBC{4}(x4,y4))];

iiv = [ii{1}; ii{2}; ii{3}; ii{4}; ii{5}];

oov = S*iiv + c;

ubv = (iiv - oov) / (2i*eta);

u = zeros(size(sd.rhs));

u(sd.idx_boundary)  = ubv;
u(sd.idx_interior)  = sd.A \ (sd.rhs(sd.idx_interior) - sd.B * ubv);

uex = u_exact(sd.px, sd.py);

fprintf("div: %i\n", div);
fprintf("Rel L2 error = %.15e\n", norm(u-uex)/norm(uex));

n = 2^div + 1;

figure;
surf(reshape(sd.px,n,n), reshape(sd.py,n,n), reshape(real(u),n,n));
view(3);