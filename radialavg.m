function R = radialavg(cor,N)
% find the radial average of the NxN correlation function cor

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