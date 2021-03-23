%-------------------------------------------------------------------------------
% SET PARAMETERS
%-------------------------------------------------------------------------------
% kT, rescaled temperature
kT = 2.5;
% N, linear lattice size
N = 80;
% J, coupling strength
J = 1; % (change sign for antiferromagnetic coupling!)
% numTimePoints, number of update steps (use large multiple of N^2 for Metropolis)
numTimePoints = 200*N^2;
% reInitialize, whether to generate a new initial condition (or continue from previous)
reInitialize = true;
% p, average proportion of initial +1 spins
p = 0.5; % (0.5 for random initial condition)
% samplingMethod, 'HeatBath', 'Metropolis' or 'Wolff'
samplingMethod = 'Metropolis';
% timeLag
timeLag = 0.1; % slow down plotting

%-------------------------------------------------------------------------------
% Generate a random initial configuration.
%-------------------------------------------------------------------------------
% Comment out to keep sampling with the previous configuration.
if reInitialize
    grid = sign(p-rand(N)); % random initial configuration
end

%-------------------------------------------------------------------------------
% Run the sampling algorithm
%-------------------------------------------------------------------------------
switch samplingMethod
case 'HeatBath'
    doHeatBath = true;
    grid = MetropolisSample(N,kT,J,numTimePoints,grid,doHeatBath,timeLag);
case 'Metropolis'
    doHeatBath = false;
    grid = MetropolisSample(N,kT,J,numTimePoints,grid,doHeatBath,timeLag);
case 'Wolff'
    numTimePoints = N^2; % run with t update steps
    grid = WolffSample(N,kT,J,numTimePoints,grid);
otherwise
    error('Description');
end
% (return a matrix of spin values: the final state)

%-------------------------------------------------------------------------------
% Compute final magnetization density and energy density
M = sum(sum(grid))/numel(grid);
E = IsingEnergy(grid,J);

f = figure(2);
f.Color = 'w';

% Plot correlation function
subplot(121)
cor = CorrelationFun(grid);
xlabel('x')
ylabel('y')
zlabel('Corr')
surf(cor);
title('Correlation')

% Plot Radial Average
subplot(122)
R = RadialAverage(cor,N);
plot(R);
xlabel('angle')
ylabel('Correlation')
title('Radial average')
