function newGrid = CoarseGrain(grid)
% 3 x 3 majority-vote coarse-graining of an input grid
%-------------------------------------------------------------------------------

% Check square grid
assert(size(grid,1)==size(grid,2))
N = size(grid,1);
% Check 3N x 3N grid for some integer N
assert(rem(N,3)==0)

newGrid = zeros(N/3,N/3);

for i = 1:N/3
    ii = (i-1)*3+1;
    for j = 1:N/3
        jj = (j-1)*3+1;
        % Get the local 3x3 grid:
        localGrid = grid(ii:ii+2,jj:jj+2);
        % Determine majority spin:
        if mean(localGrid > 0)
            newGrid(i,j) = 1;
        else
            newGrid(i,j) = -1;
        end
    end
end

end
