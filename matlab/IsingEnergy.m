function energy = IsingEnergy(grid,J)
% Compute the energy density of a spin configuration.

neighbors = circshift(grid,[0 1]) + circshift(grid,[0 -1]) + ...
            circshift(grid,[1 0]) + circshift(grid,[-1 0]);
% Each bond appears twice in the sum, hence the 1/2:
energy = -J*sum(sum(grid.*neighbors))/(2*numel(grid));

end
