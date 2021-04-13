function cor = CorrelationFun(A,doNorm)
% Computes the correlation function of A
% INPUTS: - A: a square (spin) grid.
%         - doNorm: (logical), whether to subtract the lattice-average M squared
% Output is centered at size(A)/2 because of periodic boundary conditions
%-------------------------------------------------------------------------------

% Set default:
if nargin < 2
    doNorm = true;
end

cor = conv2(A,repmat(rot90(A,2),2,2),'same')/numel(A);
if doNorm
    cor = cor - mean(A(:))^2;
end
cor = circshift(cor,ceil(size(A)/2));

end
