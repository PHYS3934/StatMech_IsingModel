% Simulate the 2D Ising model, H = -J * (sum over nearest-neighbor pairs <ij> of s_i*s_j),
% on an N x N lattice with periodic boundary conditions, then plot spin-spin correlations
% and the energy and magnetization over time.

%-------------------------------------------------------------------------------
% SET PARAMETERS
%-------------------------------------------------------------------------------
% kT, rescaled temperature (the critical temperature is kT_c = 2J/log(1+sqrt(2)), about 2.269 for J = 1)
kT_c = 2/log(1+sqrt(2));
kT = 0.85*kT_c;
% N, linear lattice size
N = 100;
% J, coupling strength (change sign for antiferromagnetic coupling!)
J = 1;
% reInitialize, whether to generate a new initial condition (or continue from previous)
reInitialize = true;
% p, average proportion of initial +1 spins
p = 0.5; % (0.5 for random initial condition)
% samplingMethod, 'HeatBath', 'Metropolis' or 'Wolff'
samplingMethod = 'HeatBath';
switch samplingMethod
case {'Metropolis','HeatBath'}
    % numTimePoints, number of single-spin updates (N^2 updates make one sweep of the lattice)
    numTimePoints = 200*N^2;
    % everyT, plot and store the energy/magnetization of the grid every everyT updates
    everyT = N^2;
    stepsPerTimeUnit = N^2;
    timeUnit = 'sweeps';
case 'Wolff'
    % numTimePoints, number of cluster flips
    % (starting from a random grid, clusters are small at first and grow as the grid orders)
    numTimePoints = 20*N;
    % everyT, plot and store the energy/magnetization of the grid every everyT cluster flips
    everyT = floor(N/5);
    stepsPerTimeUnit = 1;
    timeUnit = 'cluster flips';
end
% timeLag
timeLag = 0; % (s) pause after each round of updates (slows down plotting)
% saveVideo
saveVideo = false; % saves result to file (as a video).
% plotClusters, whether to animate the clusters and their size distribution instead of just the grid
% (try kT = kT_c to see clusters of all sizes, and a power-law distribution of cluster sizes)
plotClusters = false;

%-------------------------------------------------------------------------------
% Generate a random initial configuration
%-------------------------------------------------------------------------------
if reInitialize
    grid = sign(p - rand(N)); % random initial configuration
else
    grid = finalGrid; % from last time
end

%-------------------------------------------------------------------------------
% Run the sampling algorithm
%-------------------------------------------------------------------------------
[finalGrid,energies,magnetizations] = ...
            SampleGrid(grid,kT,J,numTimePoints,everyT,samplingMethod,timeLag,saveVideo,plotClusters);

% Time of each stored sample:
times = (0:length(energies)-1)*everyT/stepsPerTimeUnit;

%-------------------------------------------------------------------------------
% Plotting:
%-------------------------------------------------------------------------------

%------Plot spin-spin correlations------
f2 = figure(2);
f2.Color = 'w';

% Plot spin-spin correlation function (with zero displacement at the center)
subplot(121)
corrMatrix = CorrelationFun(finalGrid);
d = (1:N) - (floor(N/2)+1); % displacements matching the entries of corrMatrix
[X,Y] = meshgrid(d,d);
surf(X,Y,corrMatrix);
title('Correlation')
xlabel('x')
ylabel('y')
zlabel('Corr')
colormap('hot')
caxis([0,min(0.3,max(corrMatrix(:)))])

% Plot the radial average of spin-spin correlation function
subplot(122)
R = RadialAverage(corrMatrix,N);
plot(0:length(R)-1,R);
xlabel('Distance')
ylabel('Correlation')
title('Radial average')

%------Plot energy/magnetization variation across Markov chain evolution------
f3 = figure(3);
f3.Color = 'w';
subplot(221)
plot(times,energies)
xlabel(sprintf('time (%s)',timeUnit));
ylabel('energy per spin')
subplot(222)
histogram(energies)
xlabel('energy per spin');
ylabel('frequency')
subplot(223)
plot(times,magnetizations)
xlabel(sprintf('time (%s)',timeUnit));
ylabel('magnetization per spin')
subplot(224)
histogram(magnetizations)
xlabel('magnetization per spin');
ylabel('frequency')
