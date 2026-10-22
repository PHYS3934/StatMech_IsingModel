function cor = CorrelationFun(A,doNorm)
% Computes the spin-spin correlation function of a square spin grid, A, with periodic boundaries:
%   C(dx,dy) = mean over all sites (x,y) of A(x,y)*A(x+dx,y+dy)
% INPUTS: - A: a square (spin) grid.
%         - doNorm: (logical), whether to subtract the squared lattice-average magnetization,
%                   giving the connected correlation function
% The output is centered: zero displacement (dx = dy = 0) is at index (floor(N/2)+1, floor(N/2)+1).
%-------------------------------------------------------------------------------

% Set default:
if nargin < 2
    doNorm = true;
end

% The straightforward way loops over every displacement, shifting the grid each time:
%   cor(dx+1,dy+1) = mean(mean(A.*circshift(A,[-dx,-dy])));
% The Fourier transform gives the same result for all displacements at once, because the
% autocorrelation is the inverse transform of the power spectrum, |FFT(A)|^2
% (the Wiener-Khinchin theorem; the FFT assumes periodic boundaries, as our lattice has):
cor = real(ifft2(abs(fft2(A)).^2))/numel(A);

if doNorm
    cor = cor - mean(A(:))^2;
end

% Move zero displacement from the corner (1,1) to the center:
cor = fftshift(cor);

end
