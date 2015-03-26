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
global lambda
lambda = [1 0.8];

% TODO: How to take a Stefan parameter?
% Stefan parameter
global St
St = 1;

% Initial condition
f_init = @(z) -1 + 2*z;

% Boundary conditions
% Simple consistent conditions
global f_bottom
f_bottom = @(t) -1;

global f_top
f_top = @(t) 1;

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

%% Solving
clc

% Two grids T(z,t) -> T(x,t); z = z(x,t)

% Spatial auxiliary grid in transformed space
% [solid liquid] + additional crystallization node
N = [10 10];
xstep = 1 / sum(N);
xgrid = -0.5*xstep:xstep:(1 + 0.5*xstep);
delta_index = N(1) + 2;

% Temporal grid
tN = 100;
tmax = 0.02;
tmin = 0;
tstep = (tmax - tmin) / tN;
tgrid = tmin:tstep:tmax;

% Initial values
zcryst = zeros(1, tN+1);
zcryst(1) = 0.5;
vgrowth = zeros(1, tN+1);

% Solution matrix
T = zeros(tN, sum(N)+2);
T(1, :) = f_init(ztransform(xgrid, zcryst(1)));

%% Iterations
f_without = @(x, index) [x(1:index-1), x(index+1:end)];

maxiters = 50;
for t = 2:length(tgrid)
	
	eps = zeros(1, maxiters);
	
	vgrowth(t) = vgrowth(t-1);
	delta_new = vgrowth(t) * tstep;
	
	for i = 1:maxiters
		[a, b, c, f] = get_coeff(T(t-1, :), zcryst(t-1) + delta_new, tgrid(t), N, vgrowth(t), xstep, tstep);
		A = diag(c) + diag(a, -1) + diag(b, 1);
		F = f';
		
		Solution = prog(A, F);
		
		delta_old = delta_new;
		delta_new = Solution(delta_index);
		vgrowth(t) = delta_new / tstep;
		
		eps(i) = abs(delta_new - delta_old) / delta_new;
	end
	
	T(t, :) = f_without(Solution', delta_index);
	zcryst(t) = zcryst(t-1) + delta_new;
	
	% plot(eps)
	% pause()

end





















