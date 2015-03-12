function zgrid = ztransform( xgrid, N, zcryst )
% z-transform of grid in Stefan problem
% z is a grid, where the point of crystalization lie in the zcryst

% INPUT:
% :x: node of x-grid
% :N: numbers of grid nodes N = [solid, liquid]
% :zcryst: z-coordinate of the crystalization point

% OUTPUT:
% :z: values of corresponding node of z-grid

zgrid = zeros(1, sum(N)+2);
zgrid(1:N(1)+1) = 0.5 * xgrid(1:N(1)+1) * zcryst;
zgrid(N(1)+2:end) = 0.5 * xgrid(N(1)+2:end) * (1 - zcryst);

end

