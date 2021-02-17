function isingplot(grid,N,J,kT,M,E)
% Display the current state of the system
image((grid+1)*128);
title(sprintf('2D Ising model with %0.4g by %0.4g lattice',N,N));
xlabel(sprintf('J=%0.2f, kT = %0.2f, M = %0.3f, E = %0.3f',J,kT,M,E));
axis square; colormap bone; set(gca,'XTickLabel',[],'YTickLabel',[]);
drawnow;
end
