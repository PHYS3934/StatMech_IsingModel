function R = RadialAverage(cor,N)
% Compute the radial average of the NxN correlation function, cor
% (average out the angular dependence from a 2D connected correlation function).
% Zero displacement must be at the center, cor(floor(N/2)+1,floor(N/2)+1), as returned by CorrelationFun.
% R(r+1) is the average of cor over all displacements at distance r (to the nearest integer).
%-------------------------------------------------------------------------------

% Distance of each displacement (dx,dy) from zero displacement at the center:
d = (1:N) - (floor(N/2)+1);
[dx,dy] = meshgrid(d,d);
distance = sqrt(dx.^2 + dy.^2);

% Only go out to about N/2: beyond that, displacements wrap around the periodic boundaries
numR = ceil(N/2);
R = zeros(numR,1);

for r = 0:numR-1
    % Average over the ring of displacements with r - 0.5 <= distance < r + 0.5:
    ring = (r - 0.5 <= distance) & (distance < r + 0.5);
    R(r+1) = mean(cor(ring));
end

end
