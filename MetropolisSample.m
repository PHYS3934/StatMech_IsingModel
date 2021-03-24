function [grid,energyStore,M_store] = MetropolisSample(N,kT,J,numTimePoints,grid,doHeatBath,timeLag)
% Metropolis sampling for the 2D Ising model
%-------------------------------------------------------------------------------
if nargin < 6
    doHeatBath = false;
end
if nargin < 7
    timeLag = 0;
end

%-------------------------------------------------------------------------------
% Plot initial spin configuration
%-------------------------------------------------------------------------------
M_store = zeros(floor(numTimePoints/N^2),1);
energyStore = zeros(floor(numTimePoints/N^2),1);
M_store(1) = sum(grid(:))/numel(grid);
energyStore(1) = IsingEnergy(grid,J);
f1 = figure(1);
f1.Color = 'w';
h_Image = IsingPlot(grid,N,J,kT,M_store(1),energyStore(1));

%-------------------------------------------------------------------------------
% Evolve the Markov chain for a fixed number of steps
%-------------------------------------------------------------------------------

% Precompute the indicies adjacent to each spin index
adj = myNeighbors(1:N^2,N);

% Precompute a sequence of random spins (with a linear index)
spin = randi(N^2,numTimePoints,1);

for t = 1:numTimePoints
    % Index, s, of the spin to consider flipping:
    s = spin(t);
    if doHeatBath
        % Calculate the difference in energy between s up/down
        pUp = J*sum(grid(adj(s,:)));
        pDown = -pUp;
        z = exp(-pUp/kT) + exp(-pDown/kT);
        p = exp(-pUp/kT)/z;
        % Decide whether to set this spin up or down:
        if rand <= p
            grid(s) = -1;
        else
            grid(s) = 1;
        end
    else
        % Metropolis sampling:
        deltaE = 2*J*grid(s)*sum(grid(adj(s,:)));
        if deltaE < 0
            % Always flip
            grid(s) = -grid(s);
        else
            % Calculate the transition probability
            p = exp(-deltaE/kT);
            % Decide if a transition will occur
            if rand <= p
                grid(s) = -grid(s);
            end
        end
    end

    % Refresh display of current spin configuration every N^2 trials
    if mod(t,N^2)==0
        % Sum up our variables of interest and plot:
        M = sum(grid(:))/numel(grid);
        E = IsingEnergy(grid,J);
        h_Image = IsingPlot(grid,N,J,kT,M,E,h_Image);
        % Store for later:
        energyStore(t/N^2) = E;
        M_store(t/N^2) = M;
        % Pause if required:
        if timeLag > 0
            pause(timeLag)
        end
    end
end

end
