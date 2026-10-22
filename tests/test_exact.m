% Check the MATLAB functions against exact results on small periodic lattices.
% Run from the repo root with: run('tests/test_exact.m')
%-------------------------------------------------------------------------------

addpath(fullfile(fileparts(mfilename('fullpath')),'..','matlab'));
rng(0);

%-------------------------------------------------------------------------------
% Energy of the all-up grid: 2 bonds per spin, so the energy per spin is -2J
%-------------------------------------------------------------------------------
assert(IsingEnergy(ones(4),1) == -2);

%-------------------------------------------------------------------------------
% Correlation function and radial average against a direct calculation
%-------------------------------------------------------------------------------
for n = [8,9] % even and odd lattice sizes
    A = sign(0.5 - rand(n));
    cor = CorrelationFun(A,false);
    c = floor(n/2) + 1; % index of zero displacement
    for dx = -floor(n/2):ceil(n/2)-1
        for dy = -floor(n/2):ceil(n/2)-1
            direct = mean(mean(A.*circshift(A,[-dx,-dy])));
            assert(abs(cor(c+dx,c+dy) - direct) < 1e-12);
        end
    end
    R = RadialAverage(cor,n);
    assert(abs(R(1) - 1) < 1e-12); % distance 0: <s^2> = 1
end
fprintf(1,'CorrelationFun and RadialAverage match the direct calculation\n');

%-------------------------------------------------------------------------------
% Clusters connect nearest neighbors only, so diagonal spins are separate clusters
% (ClusterSizeStats needs the Image Processing Toolbox)
%-------------------------------------------------------------------------------
if license('test','Image_Toolbox') && ~isempty(ver('images'))
    f = figure('Visible','off');
    assert(isequal(sort(ClusterSizeStats([1,-1;-1,1])),[1,1,1,1]));
    assert(isequal(sort(ClusterSizeStats([1,1,-1;1,-1,-1;-1,-1,-1])),[3,6]));
    close(f);
    fprintf(1,'ClusterSizeStats finds the expected clusters\n');
else
    fprintf(1,'Skipping ClusterSizeStats (no Image Processing Toolbox)\n');
end

%-------------------------------------------------------------------------------
% Wolff sampling on a 4x4 lattice at kT = 2.5 against exact enumeration
% (<|m|> = 0.7647, energy per spin = -1.3791)
%-------------------------------------------------------------------------------
N = 4; J = 1; kT = 2.5;
p = 1 - exp(-2*J/kT);
adj = myNeighbors(1:N^2,N);
grid = ones(N);
T = 30000;
absM = zeros(T,1);
E = zeros(T,1);
for t = 1:T
    cluster = WolffIteration(N,p,grid,adj);
    grid(cluster) = -grid(cluster);
    absM(t) = abs(mean(grid(:)));
    E(t) = IsingEnergy(grid,J);
end
fprintf(1,'Wolff <|m|> = %.4f (exact 0.7647), energy = %.4f (exact -1.3791)\n',...
            mean(absM(1001:end)),mean(E(1001:end)));
% Tolerances are about 4 run-to-run standard deviations
assert(abs(mean(absM(1001:end)) - 0.7647) < 0.012);
assert(abs(mean(E(1001:end)) + 1.3791) < 0.03);

fprintf(1,'All tests passed\n');
