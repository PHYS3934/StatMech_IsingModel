function adj = myNeighbors(s,N)
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
