%% Test: compute tridiagonal matrix algorithm coefficients
% 9 nodes
clc

N = [4 4];
xstep = 1 / sum(N);
xgrid = 0:xstep:1;

T = zeros(2, sum(N)+2);

zcryst = 0.5;
vgrowth = 0;
zgrid = ztransform(xgrid, zcryst);

T(1, :) = f_init(ztransform(-0.5*xstep:xstep:(1 + 0.5*xstep), zcryst));
time = 1;

[a, b, c, f] = get_coeff(T(1, :), zcryst, time, N, vgrowth, xstep, time);

% Matrix assembly
A = diag(c) + diag(a, -1) + diag(b, 1);
F = f';