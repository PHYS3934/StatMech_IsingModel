N = 20; % grid size
grid = sign(0.5-rand(N)); % random initial configuration
f = figure('color','w');
h_image = imagesc(grid);
colormap('bone');
cB = colorbar();
caxis([-1,3])
cB.Limits = [-1,3];
cB.Ticks = -1:1:3;
cB.TickLabels = {'down','','up','seed','cluster'};
axis('square')

%-------------------------------------------------------------------------------
% Compute cluster growth probability
J = 1; kT = 3;
p = 1 - exp(-2*J/kT);
fprintf(1,'Iteratively add neighbor to cluster with probability %.3f\n',p);

%-------------------------------------------------------------------------------
adj = myNeighbors(1:N^2,N);
reply = '';
while isempty(reply)
    [C,seed] = WolffIteration(N,p,grid,adj);
    theNewSpin = -grid(seed);
    grid(seed) = 2;
    h_image.CData = grid;
    reply = input('What do you think of my seed?','s');
    grid(C) = 3;
    h_image.CData = grid;
    reply = input('What do you think of my cluster?','s');
    grid(C) = theNewSpin;
end
