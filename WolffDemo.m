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
ylabel('Corr0')
ylim([-0.1,1])
title(r)

%-------------------------------------------------------------------------------
% Compute cluster growth probability (based on J/kT)
p = 1 - exp(-2*J/kT);
fprintf(1,'Iteratively add neighbor to cluster with probability %.3f\n',p);

%-------------------------------------------------------------------------------
adj = myNeighbors(1:N^2,N);
reply = '';
subplot(5,1,1:4)
while isempty(reply)
    [C,seed] = WolffIteration(N,p,grid,adj);
    r = [r; corr(grid(:),grid0(:))];
    h_r.XData = 1:length(r);
    h_r.YData = r;

    theNewSpin = -grid(seed);
    grid(seed) = 2;
    h_image.CData = grid;
    title('SEED')
    reply = input('What do you think of my seed?','s');
    grid(C) = 3;
    grid(seed) = 2;
    h_image.CData = grid;
    title('CLUSTER')
    reply = input('What do you think of my cluster?','s');
    grid(C) = theNewSpin;

end
