function h_Image = IsingPlot(grid,N,J,kT,M,E,h_Image)
% Display the current state of the system
%-------------------------------------------------------------------------------

% Set default (set h_Image empty to initialize)
if nargin < 6
    N = [];
end
if nargin < 7
    h_Image = [];
end

%-------------------------------------------------------------------------------
if isempty(h_Image)
    % Initialize:
    h_Image = imagesc(grid);
    title(sprintf('2D Ising model with %0.4g by %0.4g lattice',N,N));
    axis('square');
    colormap([0,0,0;1,1,1]);
    ax = gca();
    ax.XTick = 1:length(grid);
    ax.YTick = 1:length(grid);
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    caxis([-1,1])
else
    % Update
    h_Image.CData = grid;
end

xlabel(sprintf('J = %0.2f, kT = %0.2f, M = %0.3f, E = %0.3f',...
                J,kT,M,E));
drawnow();

end
