function energy = IsingEnergy(grid,J)
% Compute the energy per spin of a spin configuration,
%   E/N^2 = -(J/N^2) * (sum over nearest-neighbor pairs <ij> of s_i*s_j),
% on an N x N lattice with periodic boundary conditions.
%-------------------------------------------------------------------------------

% Pair each spin with its right and down neighbors only, so that each bond is counted once
% (circshift wraps around the edges, which gives periodic boundary conditions)
right = circshift(grid,[0 -1]);
down = circshift(grid,[-1 0]);
bondSum = sum(sum(grid.*right)) + sum(sum(grid.*down));

% Energy per spin:
energy = -J*bondSum/numel(grid);

end
