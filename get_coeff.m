function [a, b, c, f] = get_coeff(T, zcryst, time, N, vgrowth, xstep, tstep)
% Evaluation of coefficients for tridiagonal matrix algorithm

% TODO: make all parameters global (or some of them) 

% INPUT:
% :T:
% :zcryst:
% :time:
% :N:
% :vgrouth:
% :lambda:
% :xstep:
% :tstep:

% OUTPUT:
% :a: a-coefficients
% :b: ...
% :c: ...
% :f: right side coefficients

% GLOBAL:
global lambda
global St
global f_bottom
global f_top

% alpha = [alpha_solid alpha_liquid]
alpha = tstep * lambda / (2 * xstep^2);
delta = tstep / (4 * xstep);

% $\dfrac{dz}{dt}$ computation
xgrid = 0:xstep:1;
zgrid = ztransform(xgrid, zcryst);

zcryst_next = zcryst + tstep * vgrowth;
zgrid_next = ztransform(xgrid, zcryst_next);

zdiff_half = (zgrid_next - zgrid) / tstep;

% $\psi = \dfrac{dz}{dx}$ computation

% jacobians = [solid liquid]
jacobian = [2 * zcryst, 2 * (1 - zcryst)];
jacobian_next= [2 * zcryst_next, 2 * (1 - zcryst_next)];
jacobian_half = 0.5 * (jacobian + jacobian_next);

% Memory allocation
a = zeros(1, sum(N));
b = zeros(1, sum(N));
c = zeros(1, sum(N));
f = zeros(1, sum(N));

% Coefficients for inner nodes 

% TODO: use repmat
% TODO: use vector selecters
a(1:N(1)-1) = alpha(1) / jacobian_half(1) - delta * zdiff_half(1:N(1)-1);
a(N(1)+2:end) = alpha(2) / jacobian_half(2) - delta * zdiff_half(N(1)+2:end-1);
		  
b(1:N(1)-1) = alpha(1) / jacobian_half(1) + delta * zdiff_half(2:N(1));
b(N(1)+2:end) = alpha(2) / jacobian_half(2) + delta * zdiff_half(N(1)+3:end); 

c(1:N(1)-1) = jacobian_next(1) + 2 * alpha(1) / jacobian_half(1) ... 
	+ delta * (zdiff_half(1:N(1)-1) - zdiff_half(2:N(1)));
c(N(1)+2:end) = jacobian_next(2) + 2 * alpha(2) / jacobian_half(2) ...
	+ delta * (zdiff_half(N(1)+2:end-1) - zdiff_half(N(1)+3:end));

f(1:N(1)-1) = jacobian(1) * T(2:N(1)) ...
	+ (alpha(1) / jacobian_half(1)) * (T(3:N(1)+1) - 2*T(2:N(1)) + T(1:N(1)-1)) ...
	+ delta * ((zdiff_half(2:N(1)) .* (T(3:N(1)+1) + T(2:N(1)))) ... 
				- (zdiff_half(1:N(1)-1) .* (T(2:N(1)) + T(1:N(1)-1))));
            
f(N(1)+2:end) = jacobian(2) * T(N(1)+3:end-1) ...
	+ (alpha(2) / jacobian_half(2)) * (T(N(1)+4:end) - 2*T(N(1)+3:end-1) + T(N(1)+2:end-2)) ...
	+ delta * ((zdiff_half(N(1)+3:end) .* (T(N(1)+4:end) + T(N(1)+3:end-1))) ... 
				 -(zdiff_half(N(1)+2:end-1) .* (T(N(1)+3:end-1) + T(N(1)+2:end-2))));

% Coefficients for solid border node

a(N(1)) = alpha(1) / jacobian_half(1) - delta * zdiff_half(N(1));
b(N(1)) = 0;
c(N(1)) = jacobian_next(1) + (3*alpha(1) / jacobian_half(1)) + delta * zdiff_half(N(1));
f(N(1)) = jacobian(1) * T(N(1)+1) - (2*alpha(1) / jacobian_half(1)) * T(N(1)+1) ...
	- (alpha(1) / jacobian_half(1)) * (T(N(1)+1) - T(N(1))) ...
	- delta * zdiff_half(N(1)+1) * (T(N(1)+1) + T(N(1)));

% Coefficients for liquid border node

a(N(1)+1) = 0;
b(N(1)+1) = alpha(2) / jacobian_half(2) + delta * zdiff_half(N(1)+2);
c(N(1)+1) = jacobian_next(2) + (3*alpha(2) / jacobian_half(2)) - delta * zdiff_half(N(1)+2);
f(N(1)+1) = jacobian(2) * T(N(1)+2) - (2*alpha(2) / jacobian_half(2)) * T(N(1)+2) ...
	+ (alpha(2) / jacobian_half(2)) * (T(N(1)+3) - T(N(1)+2)) ...
	+ delta * zdiff_half(N(1)+2) * (T(N(1)+3) + T(N(1)+2));

% Coefficients for crystalization node

lambda_sl = lambda(1) / lambda(2);
beta = tstep / xstep;

af = (lambda_sl * beta) / jacobian_half(1);
bf = beta / jacobian_half(2);
cf = -St;
ff = (lambda_sl * beta / jacobian_half(1)) * T(N(1)+1) ...
	+ (beta / jacobian_half(2)) * T(N(1)+2);

a = [a(1:N(1)) af a(N(1)+1:end)];
b = [b(1:N(1)) bf b(N(1)+1:end)];
c = [c(1:N(1)) cf c(N(1)+1:end)];
f = [f(1:N(1)) ff f(N(1)+1:end)];

% Coefficients for boundary nodes
a = [a 0.5];
b = [0.5 b];
c = [0.5 c 0.5];
f = [f_bottom(time) f f_top(time)];
% OR f = [f_top(time) f f_bottom(time)]; ???
f = -f(end:-1:1);

end

