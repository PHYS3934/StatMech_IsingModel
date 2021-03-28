function [grid,energyStore,M_store] = SampleGrid(N,kT,J,numTimePoints,grid,sampleHow,timeLag)
% Metropolis sampling for the 2D Ising model
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
% Precompute the indicies adjacent to each spin index
adj = myNeighbors(1:N^2,N);

switch sampleHow
case {'HeatBath','Metropolis'}
    % Precompute a sequence of random spins (with a linear index)
    spin = randi(N^2,numTimePoints,1);
    everyT = N^2;
case 'Wolff'
    p = 1-exp(-2*J/kT);
    everyT = N;
end

%-------------------------------------------------------------------------------
% Plot initial spin configuration
%-------------------------------------------------------------------------------
M_store = zeros(floor(numTimePoints/everyT),1);
energyStore = zeros(floor(numTimePoints/everyT),1);
M_store(1) = sum(grid(:))/numel(grid);
energyStore(1) = IsingEnergy(grid,J);
f1 = figure(1);
f1.Color = 'w';
h_Image = IsingPlot(grid,N,J,kT,M_store(1),energyStore(1));

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
        deltaE = 2*J*grid(s)*sum(grid(adj(s,:)));
        if deltaE < 0
            % Always flip to lower energy
            grid(s) = -grid(s);
        else
            % Calculate the transition probability
            p = exp(-deltaE/kT);
            % Decide if a transition will occur
            if rand <= p
                grid(s) = -grid(s);
            end
        end
    case 'Wolff'
        % Identify clusters and flip them t times using the Wolff algorithm
    	C = OneWolff(N,p,grid,adj);
        grid(C) = -grid(C);
    otherwise
        error('Unknown sampling method ''%s''',samplingMethod);
    end

    % Refresh display of current spin configuration every N^2 trials
    if mod(t,everyT)==0
        % Sum up our variables of interest and plot:
        M = sum(grid(:))/numel(grid);
        E = IsingEnergy(grid,J);
        h_Image = IsingPlot(grid,N,J,kT,M,E,h_Image);
        % Store for later:
        energyStore(t/everyT) = E;
        M_store(t/everyT) = M;
        % Pause if required:
        if timeLag > 0
            pause(timeLag)
        end
    end
end

%-------------------------------------------------------------------------------
% Wolff iteraction
%-------------------------------------------------------------------------------
function C = OneWolff(N,p,grid,adj)
	% Find a cluster, C, according the the Wolff sampling rule

	i = randi(N^2); % random seed spin
	C = i; % the cluster
	F = i; % the frontier of spins
	s = grid(i); % seed spin direction
	Ci = zeros(N^2,1); % indicator function for cluster elements

    while ~isempty(F)
        F = adj(F,:); % Compute the new neighboring spins
        F = F(grid(F(:)) == s);% only choose ones parallel to the seed spin
        % find elements that aren't in the cluster
        Fi = zeros(N^2,1);% indicator function for the frontier spins
        Ci(C) = 1;
        Fi(F) = 1;
        F = find(Fi - Ci > 0);
        F = F(rand(1,length(F)) < p); % keep spins only with probability p
        C(end+1:end+length(F)) = F; % add to cluster
	end

end

end
