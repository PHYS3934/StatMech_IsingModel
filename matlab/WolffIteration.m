function [C,i] = WolffIteration(N,p,grid,adj)
% Find a cluster, C, according the the Wolff sampling rule
%-------------------------------------------------------------------------------

i = randi(N^2); % random seed spin
C = i; % the cluster
F = i; % the frontier of spins
s = grid(i); % seed spin direction
Ci = zeros(N^2,1); % indicator function for cluster elements

while ~isempty(F)
    F = adj(F,:); % Compute the new neighboring spins
    F = F(grid(F(:)) == s); % only choose ones parallel to the seed spin
    % Find elements that aren't in the cluster
    Fi = zeros(N^2,1); % indicator function for the frontier spins
    Ci(C) = 1;
    Fi(F) = 1;
    F = find(Fi - Ci > 0);
    F = F(rand(1,length(F)) < p); % keep spins only with probability p
    C(end+1:end+length(F)) = F; % add to cluster
end

end
