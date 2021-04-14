function [grid,energyStore,M_store] = SampleGrid(grid,kT,J,numTimePoints,everyT,sampleHow,timeLag)
% Sampling algorithms for the 2D Ising model
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 6
    sampleHow = 'Metropolis';
end
if nargin < 7
    timeLag = 0;
end

%-------------------------------------------------------------------------------
% Preparation
%-------------------------------------------------------------------------------
N = size(grid,1);

% Precompute the indicies adjacent to each spin index
adj = myNeighbors(1:N^2,N);

switch sampleHow
case {'HeatBath','Metropolis'}
    % Precompute a sequence of random spins (with a linear index)
    spin = randi(N^2,numTimePoints,1);
case 'Wolff'
    p = 1 - exp(-2*J/kT);
end

%-------------------------------------------------------------------------------
% Plot the initial spin configuration
%-------------------------------------------------------------------------------
% M_store = zeros(floor(numTimePoints/everyT),1);
% energyStore = zeros(floor(numTimePoints/everyT),1);
% M_store(1) = mean(grid(:));
% energyStore(1) = IsingEnergy(grid,J);
f1 = figure(1);
f1.Color = 'w';
% h_Image = IsingPlot(grid,N,J,kT,M_store(1),energyStore(1));

%-------------------------------------------------------------------------------
% Evolve the Markov chain for numTimePoints iterations
%-------------------------------------------------------------------------------
for t = 1:numTimePoints
    switch sampleHow
    case 'HeatBath'
        % Index, s, of the spin to consider flipping:
        s = spin(t);
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

    case 'Metropolis'
        % Index, s, of the spin to consider flipping:
        s = spin(t);
        % Compute the change in energy from flipping this spin:
        deltaE = 2*J*grid(s)*sum(grid(adj(s,:)));
        if deltaE < 0
            % Always flip to lower energy
            grid(s) = -grid(s);
        else
            % Calculate the transition probability
            p = exp(-deltaE/kT);
            % A transition to higher energy occurs with probability p:
            if rand <= p
                grid(s) = -grid(s);
            end
        end

    case 'Wolff'
        % Identify a cluster to flip using the Wolff algorithm
    	C = WolffIteration(N,p,grid,adj);
        grid(C) = -grid(C);

    otherwise
        error('Unknown sampling method ''%s''',samplingMethod);
    end

    % Refresh display of current spin configuration every N^2 trials
    if mod(t,everyT)==0
        ClusterSizeStats(grid);
        drawnow()
        % Sum up our variables of interest and plot:
        % M = sum(grid(:))/numel(grid);
        % E = IsingEnergy(grid,J);
        % h_Image = IsingPlot(grid,N,J,kT,M,E,h_Image);
        % Store for later:
        % energyStore(t/everyT) = E;
        % M_store(t/everyT) = M;
        % Pause if required:
        if timeLag > 0
            pause(timeLag)
        end
    end
end


end
