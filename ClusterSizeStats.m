function clusterSizes = ClusterSizeStats(grid)
% Get statistics on cluster sizes
%-------------------------------------------------------------------------------

CC1 = bwconncomp(grid > 0);
SS1 = regionprops(CC1,'Area');
CC2 = bwconncomp(grid < 0);
SS2 = regionprops(CC2,'Area');
clusterSizes = [[SS1.Area],[SS2.Area]]; % cluster areas

%-------------------------------------------------------------------------------
% Logarithmic binning
%-------------------------------------------------------------------------------
numBins = ceil(min(length(unique(clusterSizes)/5),15));

% log10-spaced bin edges:
binEdges = logspace(log10(min(clusterSizes*0.9999)),log10(max(clusterSizes*1.0001)),numBins);

% Bin the data using custom bin edges:
[N,binEdges] = histcounts(clusterSizes,binEdges);

% Bin centers as middle points between bin edges:
binCenters = mean([binEdges(1:end-1);binEdges(2:end)]);

% Convert counts to probabilities:
Nnorm = N/sum(N);

%-------------------------------------------------------------------------------
% PLOTTING
%-------------------------------------------------------------------------------
% f = figure('color','w');
subplot(1,3,1)
imagesc(grid)
axis('square')
colormap([0,0,0;1,1,1])
subplot(1,3,2)
L1 = labelmatrix(CC1);
L2 = labelmatrix(CC2);
RGB1 = label2rgb(L1,'jet','k','noshuffle');
RGB2 = label2rgb(L2,'lines','k','noshuffle');
imagesc(RGB1+RGB2)
axis('square')
subplot(1,3,3)
loglog(binCenters,Nnorm,'o-k')
xlabel('Cluster size')
ylabel('Probability')

end
