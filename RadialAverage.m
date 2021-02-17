function R = RadialAverage(cor,N)
% find the radial average of the NxN correlation function cor
% (average out the angular dependence from a 2D connected correlation function)

L = ceil(N/2);
c = (1:N)-L;
[X,Y] = meshgrid(c,c);
[~,rho] = cart2pol(X,Y);
rbins = -.5:1:L; % only go to N/2, because of periodic boundaries
R = zeros(L,1);
for j = 1:L,
    r = rbins(j);
    ring = (r <= rho).*(rho < r+1);
    R(j) = sum(sum(cor.*ring))/nnz(ring);
end

R(isnan(R)) = 0;

end
