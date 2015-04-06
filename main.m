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
St = 0.5;

% Initial condition
% f_init = @(z) -1 + 2*z;
f_init = @(z) 0;

% Boundary conditions
% Simple consistent conditions
global f_bottom
f_bottom = @(t) -1;

global f_top
f_top = @(t) 1;

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
tN = 50;
tmax = 1;
tmin = 0;
tstep = (tmax - tmin) / tN;
tgrid = tmin:tstep:tmax;

% Initial values
zcryst = zeros(1, tN+1);
zcryst(1) = 0.2;
vgrowth = zeros(1, tN+1);
vgrowth(1) = 0;

%% Iterations
f_without = @(x, index) [x(1:index-1), x(index+1:end)];

% Solution matrix
T = zeros(tN, sum(N)+2);
T(1, :) = f_init(ztransform(xgrid, zcryst(1)));

maxiters = 50;
for t = 1:length(tgrid)
	
	eps = zeros(1, maxiters);
	delta_new = vgrowth(t) * tstep;
	
	for i = 1:maxiters
		[a, b, c, f] = get_coeff(T(t, :), zcryst(t) + delta_new, tgrid(t), N, vgrowth(t), xstep, tstep);
		
		A = diag(c) + diag(a, -1) + diag(b, 1);
		F = f';
		
		Solution = prog(A, F);
		
		delta_old = delta_new;
		delta_new = Solution(delta_index);
		vgrowth(t) = delta_new / tstep;
		
		T(t, :) = f_without(Solution', delta_index);
		
		% eps(i) = abs(delta_new - delta_old) / delta_new;
        eps(i) = abs(delta_new - delta_old);
	end
	
	T(t+1, :) = T(t, :);
	zcryst(t+1) = zcryst(t) + delta_new;
	vgrowth(t+1) = vgrowth(t);
	% plot(eps)
	% pause()

end





















