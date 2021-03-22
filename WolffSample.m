function grid = WolffSample(N,kT,J,t,grid)
% Wolff cluster flip sampling for the 2D Ising model

% plot initial configuration
M = sum(sum(grid))/N^2;
E = IsingEnergy(grid,J);

figure(1);
IsingPlot(grid,N,J,kT,M,E);

% precompute the indicies adjacent to each spin index
adj = myNeighbors(1:N^2,N);
p = 1-exp(-2*J/kT);

for j = 1:t,
    % Identify clusters and flip them t times using the Wolff algorithm
	C = OneWolff(N,p,grid,adj);
	grid(C) = -1*grid(C);
    % Plot the relevant variables
    if mod(j,N^2)==0,
        M = sum(sum(grid))/N^2;
        E = IsingEnergy(grid,J);
        IsingPlot(grid,N,J,kT,M,E);
    end
end

end

function C = OneWolff(N,p,grid,adj)
	% Find a cluster, C, according the the Wolff sampling rule

	i = randi(N^2); % random seed spin
	C = i; % the cluster
	F = i; % the frontier of spins
	s = grid(i); % seed spin direction
	Ci = zeros(N^2,1); % indicator function for cluster elements

	while ~isempty(F),
	   F = adj(F,:); % Compute the new neighboring spins
	   F = F(grid(F(:)) == s);% only choose ones parallel to the seed spin
	   % find elements that aren't in the cluster
	   Fi = zeros(N^2,1);% indicator function for the frontier spins
	   Ci(C) = 1; Fi(F)=1;
	   F = find(Fi-Ci > 0);
	   F = F(rand(1,length(F)) < p); % keep spins only with probability p
	   C(end+1:end+length(F)) = F;% add to cluster
	end

end
