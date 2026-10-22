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
    F = F(:); % keep duplicates: each entry is one bond, with its own trial
    Ci(C) = 1;
    F = F(grid(F) == s & Ci(F) == 0); % only parallel spins not already in the cluster
    % Activate each bond with probability p, then remove duplicate sites
    % (a site touching k cluster spins joins with probability 1 - (1-p)^k)
    F = unique(F(rand(length(F),1) < p));
    C(end+1:end+length(F)) = F; % add to cluster
end

end
