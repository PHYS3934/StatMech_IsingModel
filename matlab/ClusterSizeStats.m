function clusterSizes = ClusterSizeStats(grid)
% Get statistics on cluster sizes, where a cluster is a connected region of parallel spins
% (nearest neighbors only; clusters are not joined across the periodic boundaries).
% Plots the grid, its clusters (colored by direction and size), and the distribution of
% cluster sizes; returns the cluster sizes.
% Call it repeatedly as a grid evolves to animate the clusters (see the plotClusters option of SampleGrid).
% Needs the Image Processing Toolbox (for bwconncomp).
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Find the clusters
%-------------------------------------------------------------------------------
clusterSizes = [];
sizeMap = zeros(size(grid)); % the size of the cluster that each spin belongs to
for spin = [1,-1]
    % Connected regions of this spin direction (4: connect nearest neighbors only, not diagonals):
    CC = bwconncomp(grid == spin,4);
    sizes = cellfun(@numel,CC.PixelIdxList);
    % Give each spin the size of its cluster:
    for k = 1:CC.NumObjects
        sizeMap(CC.PixelIdxList{k}) = sizes(k);
    end
    clusterSizes = [clusterSizes,sizes];
end

%-------------------------------------------------------------------------------
% Logarithmic binning
%-------------------------------------------------------------------------------
% log10-spaced bin edges from 1 to the number of spins (the largest possible cluster).
% Using the same bins for every grid of a given size lets us compare grids directly.
numBins = 20;
binEdges = logspace(0,log10(numel(grid)+1),numBins+1);

% Count the clusters in each bin, and convert counts to probabilities:
counts = histcounts(clusterSizes,binEdges);
probabilities = counts/sum(counts);

% Bin centers, as the geometric mean of the bin edges (the midpoint on a log scale):
binCenters = sqrt(binEdges(1:end-1).*binEdges(2:end));

%-------------------------------------------------------------------------------
% Color each spin by its cluster's direction and size:
% up clusters from red (small) to yellow (large), down clusters from blue (small) to green (large).
% Neighboring clusters always point in opposite directions, so they always get different colors.
%-------------------------------------------------------------------------------
logSize = log(sizeMap)/log(numel(grid)); % from 0 (one spin) to 1 (the whole grid)
colorsUp = ind2rgb(round(255*logSize) + 1,autumn(256));
colorsDown = ind2rgb(round(255*logSize) + 1,winter(256));
isUp = repmat(grid > 0,1,1,3);
clusterColors = colorsDown;
clusterColors(isUp) = colorsUp(isUp);

%-------------------------------------------------------------------------------
% PLOTTING
%-------------------------------------------------------------------------------
% The spin configuration:
subplot(1,3,1)
imagesc(grid)
colormap(gca,[0,0,0;1,1,1])
caxis([-1,1])
axis('square')
set(gca,'XTick',[],'YTick',[])
title('Spins')

% Each cluster colored by its direction and size:
subplot(1,3,2)
image(clusterColors)
axis('square')
set(gca,'XTick',[],'YTick',[])
title('Clusters, colored by size')
xlabel({'up: red (small) to yellow (large)','down: blue (small) to green (large)'})

% The distribution of cluster sizes (empty bins are left out of the log-log plot):
subplot(1,3,3)
show = probabilities > 0;
loglog(binCenters(show),probabilities(show),'o-k')
xlim([0.7,numel(grid)*1.5])
axis('square')
xlabel('Cluster size')
ylabel('Probability')
title('Cluster-size distribution')

end
