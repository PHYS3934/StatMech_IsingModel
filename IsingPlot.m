function h_Image = IsingPlot(grid,N,J,kT,M,E,h_Image)
% Display the current state of the system
%-------------------------------------------------------------------------------

% Set default (set h_Image empty to initialize)
if nargin < 7
    h_Image = [];
end

%-------------------------------------------------------------------------------
% Update:
if isempty(h_Image)
    % Initialize:
    h_Image = imagesc((grid+1)*128);
    title(sprintf('2D Ising model with %0.4g by %0.4g lattice',N,N));
    axis('square');
    colormap('bone');
    ax = gca();
    ax.XTickLabel = [];
    ax.YTickLabel = [];
else
    % Update
    h_Image.CData = (grid+1)*128;
end

xlabel(sprintf('J = %0.2f, kT = %0.2f, M = %0.3f, E = %0.3f',J,kT,M,E));
drawnow;

end
