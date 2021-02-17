function energy = isingenergy(grid,J)

neighbors = circshift(grid,[0 1]) + circshift(grid,[0 -1]) + ... 
            circshift(grid,[1 0]) + circshift(grid,[-1 0]);
energy = -J*sum(sum(grid.*neighbors))/numel(grid);

end