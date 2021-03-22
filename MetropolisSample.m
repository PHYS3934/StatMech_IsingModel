function grid = MetropolisSample(N,kT,J,t,grid,timeLag)
% Metropolis sampling for the 2D Ising model
%-------------------------------------------------------------------------------
if nargin < 6
    timeLag = 0;
end

% Plot initial spin configuration
M = sum(sum(grid))/numel(grid);
E = IsingEnergy(grid,J);
figure(1);
IsingPlot(grid,N,J,kT,M,E);

% precompute the indicies adjacent to each spin index
adj = myNeighbors(1:N^2,N);

% Pick a sequence of random spins (with a linear index)
spin = randi(N^2,t,1);

% Evolve the system for a fixed number of steps
for i = 1:t,
    % location of ith spin
    s = spin(i);
    % Calculate the change in energy from flipping s
    deltaE = 2*J*grid(s)*sum(grid(adj(s,:)));
    % Calculate the transition probability
    p = exp(-deltaE/kT);
    % Decide if a transition will occur
    if rand <= p
        grid(s) = -1*grid(s);
    end
    % Refresh display of current spin configuration every N^2 trials
    if mod(i,N^2)==0,
        % Sum up our variables of interest and plot
        M = sum(sum(grid))/numel(grid);
        E = IsingEnergy(grid,J);
        IsingPlot(grid,N,J,kT,M,E);
        if timeLag > 0
            pause(timeLag)
        end
    end
end

end
