function grid = metropolis(N,kT,J,t,grid)
% Metropolis sampling for the 2D Ising model

% plot initial configuration
M = sum(sum(grid))/numel(grid);
E = isingenergy(grid,J);
figure(1);
isingplot(grid,N,J,kT,M,E);

% precompute the indicies adjacent to each spin index
adj = neighbors(1:N^2,N);

% Pick a sequence of random spins (with a linear index)
spin = randi(N^2,t,1);

% Evolve the system for a fixed number of steps
for i=1:t,
    % location of ith spin
    s = spin(i);
    % Calculate the change in energy of flipping s
    dE = 2*J*grid(s)*sum(grid(adj(s,:)));
    % Calculate the transition probability
    p = exp(-dE/kT);
    % Decide if a transition will occur
    if rand <= p, grid(s) = -1*grid(s); end
    % Refresh display of current spin configuration every N^2 trials
    if mod(i,N^2)==0,
        % Sum up our variables of interest and plot
        M = sum(sum(grid))/numel(grid);
        E = isingenergy(grid,J);
        isingplot(grid,N,J,kT,M,E);
    end
end

end

function adj = neighbors(s,N)
% take a list of linear indices s and return the linear indices of the
% neighbors of s on an N by N grid with periodic boundary conditions.

s = s-1; % index by zero
adj = zeros(length(s),4);

% s = r*N+c;
r = floor(s/N);
c = rem(s,N);

adj(:,1) = mod(r+1,N)*N+c;  %down
adj(:,2) = mod(r-1,N)*N+c;  %up
adj(:,3) = r*N+mod(c+1,N);  %right
adj(:,4) = r*N+mod(c-1,N);  %left

adj = adj+1; % index by one again

end