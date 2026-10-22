function grid = ProppWilson(N,kT,J)
% Exact sampling of the Ising Model via the Propp-Wilson algorithm
%-------------------------------------------------------------------------------

% Initialization: chain 1 starts all up, chain 2 all down.
% Both chains use the same random numbers, so once they agree, they agree forever (coalescence).
grid = ones(N,N,2); grid(:,:,2) = -grid(:,:,2);
scale = 2*J/kT;
p = 1./(1+exp(-scale*[-4:2:4])); % precompute allowed probability values
k = 1;
V = randi(N,2,2^(k-1));
U = rand(1);
while (~isequal(grid(:,:,1),grid(:,:,2)))
    fprintf(1,'k = %u\n',k); % so you know how long it's running
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
    % Plot both chains, each with its own magnetization and energy per spin:
    for c = 1:2
        subplot(1,2,c)
        M = mean(mean(grid(:,:,c)));
        E = IsingEnergy(grid(:,:,c),J);
        IsingPlot(grid(:,:,c),N,J,kT,M,E);
    end
end

grid = grid(:,:,1);

end
