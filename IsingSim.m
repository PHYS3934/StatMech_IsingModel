%-------------------------------------------------------------------------------
% SET PARAMETERS
%-------------------------------------------------------------------------------
% kT, rescaled temperature
% kT = 1;
kT = 2/log(1+sqrt(2))*0.85;
% N, linear lattice size
N = 100;
% J, coupling strength (change sign for antiferromagnetic coupling!)
J = 1;
% reInitialize, whether to generate a new initial condition (or continue from previous)
reInitialize = true;
% p, average proportion of initial +1 spins
p = 0.5; % (0.5 for random initial condition)
% samplingMethod, 'HeatBath', 'Metropolis' or 'Wolff'
samplingMethod = 'Metropolis';
switch samplingMethod
case {'Metropolis','HeatBath'}
    % numTimePoints, number of update steps (use large multiple of N^2 for Metropolis)
    numTimePoints = 200*N^2;
    % everyT, plot and store the energy/magnetization of the grid everyT iterations
    everyT = N^2;
case 'Wolff'
    % numTimePoints, number of update steps (use large multiple of N^2 for Metropolis)
    numTimePoints = 100*N;
    % everyT, plot and store the energy/magnetization of the grid everyT iterations
    everyT = 0.5*N;
end
% timeLag
timeLag = 0; % (s) pause after each round of updates (slows down plotting)
% saveVideo
saveVideo = false; % saves result to file (as a video).

%-------------------------------------------------------------------------------
% Generate a random initial configuration
%-------------------------------------------------------------------------------
% Comment out to keep sampling with the previous configuration.
if reInitialize
    grid = sign(p-rand(N)); % random initial configuration
else
    grid = finalGrid; % from last time
end

%-------------------------------------------------------------------------------
% Run the sampling algorithm
%-------------------------------------------------------------------------------
[finalGrid,energies,magnetizations] = ...
            SampleGrid(grid,kT,J,numTimePoints,everyT,samplingMethod,timeLag,saveVideo);

%-------------------------------------------------------------------------------
% Plotting:
%-------------------------------------------------------------------------------

%------Plot spin-spin correlations------
f2 = figure(2);
f2.Color = 'w';

% Plot spin-spin correlation function
subplot(121)
corrMatrix = CorrelationFun(finalGrid);
[X,Y] = meshgrid(-N/2:(N/2-1),-N/2:(N/2-1));
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
plot(R);
xlabel('Distance')
ylabel('Correlation')
title('Radial average')

%------Plot energy/magnetization variation across Markov chain evolution------
if exist('energies','var')
    f3 = figure(3);
    f3.Color = 'w';
    subplot(221)
    plot(energies)
    xlabel(sprintf('time (iterations/%u)',everyT));
    ylabel('energy')
    subplot(222)
    histogram(energies)
    xlabel('energy');
    ylabel('frequency')
    subplot(223)
    plot(magnetizations)
    xlabel(sprintf('time (iterations/%u)',everyT));
    ylabel('magnetization')
    subplot(224)
    histogram(magnetizations)
    xlabel('magnetization');
    ylabel('frequency')
end
