%% Description
% Numeriacal study of Stefan problem

% TODO description

% We suppose the spatial interval is equal [0 1].
% Spatial transformation function z(x,t) define partialy uniform
% tranformation (uniform in solid and liquid).

%% Init
clc
clear

% Thermal conduction coefficients:
% [solid liquid]
lambda = [1 0.8];

% TODO: How to take a Stefan parameter?
% Stefan parameter
St = 1;

% Initial condition
f_init = @(z) -1 + 2*z;

% Boundary conditions
% Simple consistent conditions
f_bottom = @(t) -1;
f_top = @(t) 1;

%% Test: compute tridiagonal matrix algorithm coefficients
% 9 nodes
clc

N = [4 4];
xstep = 1 / sum(N);
xgrid = 0:xstep:1;

T = zeros(2, sum(N)+2);

zcryst = 0.5;
vgrowth = 0.5;
zgrid = ztransform(xgrid, zcryst);

T(1, :) = f_init(ztransform(-0.5*xstep:xstep:(1 + 0.5*xstep), zcryst));
time = 1;

[A, B, C, F] = get_coeff(T(1, :), zcryst, time, N, vgrowth, lambda, St, xstep, time, f_bottom, f_top);

%% Solving
clc

% Two grids T(z,t) -> T(x,t); z = z(x,t)

% Spatial auxiliary grid in transformed space
% [solid liquid] + additional crystallization node
N = [50 50];
xstep = 1 / sum(N);
xgrid = (-0.5*xstep):xtep:(1 + 0.5*xstep);

% Temporal grid
tN = 100;
tmax = 1;
tstep = 1 / (tN - 1);
tgrid = 0:tstep:tmax;

% Solution matrix
T = zeros(tN, sum(N)+1);
T(1, :) = f_init(zgrid);





















