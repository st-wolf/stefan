function zgrid = ztransform( xgrid, zcryst )
% z-transform of grid in Stefan problem
% z is a grid, where the point of crystalization lie in the zcryst

% INPUT:
% :x: node of x-grid
% :N: numbers of grid nodes N = [solid, liquid]
% :zcryst: z-coordinate of the crystalization point

% OUTPUT:
% :z: values of corresponding node of z-grid

zgrid = zeros(1, length(xgrid));

zgrid(xgrid <= 0.5) = 2 * xgrid(xgrid <= 0.5) * zcryst;
zgrid(xgrid > 0.5) = 2 * xgrid(xgrid > 0.5) * (1 - zcryst) + 2*zcryst - 1;

end

