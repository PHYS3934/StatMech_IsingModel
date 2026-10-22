% Parameters:
N = 100; % grid size
J = 1;
% kT = 2.3;
kT = 2/log(1+sqrt(2))*1.0; % T_c

%-------------------------------------------------------------------------------
grid0 = sign(0.5-rand(N)); % random initial configuration
grid = grid0;
f = figure('color','w');
ax = subplot(5,1,1:4);
ax.XTick = [];
ax.YTick = [];
h_image = imagesc(grid);
cB = colorbar();
caxis([-1,3])
cB.Limits = [-1,3];
cB.Ticks = -1:1:3;
cB.TickLabels = {'down','','up','seed','cluster'};
axis('square')
colormap([1,1,1;0.5,0.5,0.5;0,0,0;[230, 85, 13]/255;[254, 230, 206]/255]);

%-------------------------------------------------------------------------------
% Correlation to original grid
subplot(5,1,5)
r = corr(grid(:),grid0(:));
h_r = plot(r,'.-k');
xlabel('Iteration')
ylabel('Corr. with initial grid')
ylim([-0.1,1])
title(r)

%-------------------------------------------------------------------------------
% Probability of adding a parallel neighbor to the cluster (activating a bond)
p = 1 - exp(-2*J/kT);
fprintf(1,'Iteratively add neighbor to cluster with probability %.3f\n',p);

%-------------------------------------------------------------------------------
% Step through the algorithm: show the seed, then its cluster, then flip the cluster
% (press Enter to continue; type anything then Enter to stop)
adj = myNeighbors(1:N^2,N);
reply = '';
subplot(5,1,1:4)
while isempty(reply)
    % Grow a cluster from a random seed spin:
    [cluster,seed] = WolffIteration(N,p,grid,adj);

    % Show the seed (grid values 2) and then the cluster (grid values 3):
    displayGrid = grid;
    displayGrid(seed) = 2;
    h_image.CData = displayGrid;
    title('SEED')
    reply = input('What do you think of my seed?','s');
    displayGrid(cluster) = 3;
    displayGrid(seed) = 2;
    h_image.CData = displayGrid;
    title('CLUSTER')
    reply = input('What do you think of my cluster?','s');

    % Flip every spin in the cluster:
    grid(cluster) = -grid(cluster);
    h_image.CData = grid;
    title('FLIPPED')

    % Track the correlation with the initial grid:
    r = [r; corr(grid(:),grid0(:))];
    h_r.XData = 1:length(r);
    h_r.YData = r;
end
