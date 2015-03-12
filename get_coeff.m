function [a, b, c, f] = get_coeff( T, zcryst, time, N, vgrowth, lambda, St, xstep, tstep, f_bottom, f_top )
% Evaluation of coefficients for tridiagonal matrix algorithm

% TODO: make all parameters global (or some of them) 

% INPUT:
% :T:
% :zcryst:
% :time:
% :N:
% :vgrouth:
% :lambda:
% :St:
% :xstep:
% :tstep:
% ...

% OUTPUT:
% :a: a-coefficients
% :b: ...
% :c: ...
% :f: right side coefficients

% alpha = [solid liquid]
alpha = tstep * lambda / (2 * xstep^2);
delta = tstep / (4 * xstep);

% $\dfrac{dz}{dt}$ computation
xgrid = (-0.5*xstep):xstep:(1 + 0.5*xstep);
zgrid = ztransform(xgrid, N, zcryst);

zcryst_next = zcryst + tstep * vgrowth;
zgrid_next = ztransform(xgrid, N, zcryst_next);

zdiff_half = (zgrid_next - zgrid) / tstep;

% $\psi = \dfrac{dz}{dx}$ computation

% jacobians = [solid liquid]
jacobian = [0.5 * zcryst 0.5 * (1 - zcryst)];
jacobian_next= [0.5 * zcryst_next 0.5 * (1 - zcryst_next)];
jacobian_half = 0.5 * (jacobian + jacobian_next);

% Memory allocation
a = zeros(1, sum(N)+2);
b = zeros(1, sum(N)+2);
c = zeros(1, sum(N)+2);
f = zeros(1, sum(N)+2);

% Coefficients for inner nodes 

% TODO: use repmat
% TODO: use vector selecters
a(2:N(1)-1) = alpha(1) / jacobian_half(1) - delta * zdiff_half(2:N(1)-1);
a(N(1)+3:sum(N)) = alpha(2) / jacobian_half(2) - delta * zdiff_half(N(1)+3:end-1);
		  
b(2:N(1)-1) = alpha(1) / jacobian_half(1) + delta * zdiff_half(3:N(1));
b(N(1)+3:sum(N)) = alpha(2) / jacobian_half(2) + delta * zdiff_half(N(1)+4:end); 

c(2:N(1)-1) = jacobian_next(1) + 2 * alpha(1) / jacobian_half(1) ... 
	+ delta * (zdiff_half(2:N(1)-1) - zdiff_half(3:N(1)));
c(N(1)+3:sum(N)) = jacobian_next(2) + 2 * alpha(2) / jacobian_half(2) ...
	+ delta * (zdiff_half(N(1)+3:end-1) - zdiff_half(N(1)+4:end));

f(2:N(1)-1) = jacobian(1) * T(2:N(1)-1) ...
	+ (alpha(1) / jacobian_half(1)) * (T(3:N(1)) - 2*T(2:N(1)-1) + T(1:N(1)-2)) ...
	+ delta * ((zdiff_half(3:N(1)) .* (T(3:N(1)) + T(2:N(1)-1))) ... 
				- (zdiff_half(2:N(1)-1) .* (T(2:N(1)-1) + T(1:N(1)-2))));
			
f(N(1)+3:sum(N)) = jacobian(2) * T(N(1)+3:sum(N)) ...
	+ (alpha(2) / jacobian_half(2)) * (T(N(1)+4:end) - 2*T(N(1)+3:end-1) + T(N(1)+2:end-2)) ...
	+ delta * ((zdiff_half(N(1)+4:end) .* (T(N(1)+4:end) + T(N(1)+3:end-1))) ... 
				 -(zdiff_half(N(1)+3:end-1) .* (T(N(1)+3:end-1) + T(N(1)+2:end-2))));

% Coefficients for solid border node

a(N(1)) = alpha(1) / jacobian_half(1) - delta * zdiff_half(N(1));
b(N(1)) = 0;
c(N(1)) = jacobian_next(1) + (3*alpha(1) / jacobian_half(1)) + delta * zdiff_half(N(1));
f(N(1)) = jacobian(1) * T(N(1)) - (4*alpha(1) / sum(jacobian_half)) * T(N(1)) ...
	- (alpha(1) / jacobian_half(1)) * (T(N(1)) - T(N(1)-1)) ...
	- delta * zdiff_half(N(1)) * (T(N(1)) + T(N(1)-1));

% Coefficients for liquid border node

a(N(1)+2) = 0;
b(N(1)+2) = alpha(2) / jacobian_half(2) + delta * zdiff_half(N(1)+2);
c(N(1)+2) = jacobian_next(2) + (3*alpha(2) / jacobian_half(2)) - delta * zdiff_half(N(1)+2);
f(N(1)+2) = jacobian(2) * T(N(1)+2) - (4*alpha(2) / sum(jacobian_half)) * T(N(1)+1) ...
	- (alpha(2) / jacobian_half(2) * (T(N(1)+2) - T(N(1)+1))) ...
	- delta * zdiff_half(N(1)+2) * (T(N(1)+2) + T(N(1)+1));

% Coefficients for crystalization node

lambda_sl = lambda(1) / lambda(2);
beta = tstep / xstep;

a(N(1)+1) = (lambda_sl * beta) / jacobian_half(1);
b(N(1)+1) = beta / jacobian_half(2);
c(N(1)+1) = -St;
f(N(1)+1) = (lambda_sl * beta / jacobian_half(1)) * T(N(1)) ...
	+ (beta / jacobian_half(2)) * T(N(1)+1);

% Coefficients for boundary nodes
a(1) = 0.5;
c(1) = 0.5;
a(end) = 0.5;
c(end) = 0.5;
f(1) = f_bottom(time);
f(end) = f_top(time);

f = -f;

end

