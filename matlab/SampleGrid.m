function [grid,energyStore,M_store] = SampleGrid(grid,kT,J,numTimePoints,everyT,sampleHow,timeLag,saveVideo)
% Evolve a 2D Ising spin configuration using a Markov chain Monte Carlo sampler,
% plotting the grid as it evolves.
%
% INPUTS: - grid: N x N matrix of spins (+1 or -1)
%         - kT: temperature (in units where Boltzmann's constant is 1)
%         - J: coupling strength
%         - numTimePoints: number of update steps. For 'HeatBath' and 'Metropolis',
%               each step updates one spin (so N^2 steps make one sweep of the lattice).
%               For 'Wolff', each step flips one cluster.
%         - everyT: plot and store the energy and magnetization every everyT steps
%         - sampleHow: 'HeatBath', 'Metropolis' (default), or 'Wolff'
%         - timeLag: pause (in seconds) after each plot (default 0)
%         - saveVideo: whether to save the animation to IsingSim.mp4 (default false)
% OUTPUTS: - grid: the final spin configuration
%          - energyStore: energy per spin at each stored time (starting with the initial grid)
%          - M_store: magnetization per spin at each stored time (starting with the initial grid)
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
if nargin < 8
    saveVideo = false;
end
if ~ismember(sampleHow,{'HeatBath','Metropolis','Wolff'})
    error('Unknown sampling method ''%s''',sampleHow);
end

%-------------------------------------------------------------------------------
% Preparation
%-------------------------------------------------------------------------------
N = size(grid,1);

% Precompute the linear indices of the four neighbors of each spin:
adj = myNeighbors(1:N^2,N);

% Probability of activating a bond in the Wolff algorithm:
pBond = 1 - exp(-2*J/kT);

if saveVideo
    v = VideoWriter('IsingSim.mp4','MPEG-4');
    open(v);
end

%-------------------------------------------------------------------------------
% Store and plot the initial spin configuration
%-------------------------------------------------------------------------------
numStored = floor(numTimePoints/everyT) + 1;
M_store = zeros(numStored,1);
energyStore = zeros(numStored,1);
M_store(1) = mean(grid(:));
energyStore(1) = IsingEnergy(grid,J);
f1 = figure(1);
f1.Color = 'w';
h_Image = IsingPlot(grid,N,J,kT,M_store(1),energyStore(1));

%-------------------------------------------------------------------------------
% Evolve the Markov chain for numTimePoints steps
%-------------------------------------------------------------------------------
for t = 1:numTimePoints
    switch sampleHow
    case 'HeatBath'
        % Pick a random spin, s (as a linear index):
        s = randi(N^2);
        % Sum of its four neighbors:
        h = sum(grid(adj(s,:)));
        % Energy of spin s pointing up or down, given its neighbors:
        E_up = -J*h;
        E_down = J*h;
        % Boltzmann probability that spin s points up:
        P_up = exp(-E_up/kT)/(exp(-E_up/kT) + exp(-E_down/kT));
        % Set spin s to up with probability P_up (regardless of its current value):
        if rand < P_up
            grid(s) = 1;
        else
            grid(s) = -1;
        end

    case 'Metropolis'
        % Pick a random spin, s (as a linear index):
        s = randi(N^2);
        % Change in energy if spin s were flipped:
        deltaE = 2*J*grid(s)*sum(grid(adj(s,:)));
        % Accept the flip with probability min(1, exp(-deltaE/kT)),
        % i.e., always flip if the energy does not increase:
        if deltaE <= 0 || rand < exp(-deltaE/kT)
            grid(s) = -grid(s);
        end

    case 'Wolff'
        % Grow a cluster using the Wolff rule, then flip all of its spins:
        cluster = WolffIteration(N,pBond,grid,adj);
        grid(cluster) = -grid(cluster);
    end

    % Store and plot the energy and magnetization every everyT steps:
    if mod(t,everyT) == 0
        k = t/everyT + 1; % index of this sample (the initial grid is sample 1)
        M_store(k) = mean(grid(:));
        energyStore(k) = IsingEnergy(grid,J);
        h_Image = IsingPlot(grid,N,J,kT,M_store(k),energyStore(k),h_Image);

        % Pause if required:
        if timeLag > 0
            pause(timeLag)
        end

        if saveVideo
            frame = getframe(f1);
            writeVideo(v,frame);
        end
    end
end

if saveVideo
    close(v);
end

end
