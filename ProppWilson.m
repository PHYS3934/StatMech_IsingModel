function grid = ProppWilson(N,kT,J)
% Exact sampling of the Ising Model via the Propp-Wilson algorithm

% Intialization
grid = ones(N,N,2); grid(:,:,2) = -grid(:,:,2);
scale = 2*J/kT;
p = 1./(1+exp(-scale*[-4:2:4])); % precompute allowed probability values
k = 1;
V = randi(N,2,2^(k-1));
U = rand(1);
while (~isequal(grid(:,:,1),grid(:,:,2)))
    display(k); % so you know how long it's running
    % Pick 2^(k-1) vertices uniformly at random
    V = [randi(N,2,2^(k-1)) V];
    % Pick 2^k uniform variables (to use for changing state)
    U = [rand(1,2^(k-1)) U];
    grid = ones(N,N,2);
    grid(:,:,2) = -grid(:,:,2);
    % Move both chains from time -2^k to time 0
    for t = 1:2^k,
        v = V(:,t)'; u = U(t);
        i = v(1); j = v(2);

        G = grid(mod(i  , N) + 1, j,:) + ...
            grid(mod(i-2, N) + 1, j,:) + ...
            grid(i, mod(j  , N) + 1,:) + ...
            grid(i, mod(j-2, N) + 1,:);

        if (u < p(G(1)/2+3))
            grid(i,j,1) = 1;
        else
            grid(i,j,1) = -1;
        end
        if (u < p(G(2)/2+3))
            grid(i,j,2) = 1;
        else
            grid(i,j,2) = -1;
        end
    end
    k = k+1;
    M = sum(sum(grid))/numel(grid(:,:,1));
    E = IsingEnergy(grid,J);
    subplot(121)
    IsingPlot(grid(:,:,1),N,J,kT,M(1),E(1));
    subplot(122)
    IsingPlot(grid(:,:,2),N,J,kT,M(1),E(1));
end

grid = grid(:,:,1);

end
