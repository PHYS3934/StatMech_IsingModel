%-------------------------------------------------------------------------------
% SET PARAMETERS
%-------------------------------------------------------------------------------
% kT, rescaled temperature
kT = 1.5; %2/log(1+sqrt(2)); %3;
% N, linear lattice size
N = 80;
% J, coupling strength (change sign for antiferromagnetic coupling!)
J = 1;
% numTimePoints, number of update steps (use large multiple of N^2 for Metropolis)
numTimePoints = 100*N^2;
% everyT, plot and store the energy/magnetization of the grid everyT iterations
everyT = N^2;
% reInitialize, whether to generate a new initial condition (or continue from previous)
reInitialize = true;
% p, average proportion of initial +1 spins
p = 0.8; % (0.5 for random initial condition)
% samplingMethod, 'HeatBath', 'Metropolis' or 'Wolff'
samplingMethod = 'Metropolis';
% timeLag
timeLag = 0; % option to slow down plotting

%-------------------------------------------------------------------------------
% Generate a random initial configuration
%-------------------------------------------------------------------------------
% Comment out to keep sampling with the previous configuration.
if reInitialize
    grid = sign(p-rand(N)); % random initial configuration
end

%-------------------------------------------------------------------------------
% Run the sampling algorithm
%-------------------------------------------------------------------------------
[finalGrid,energies,magnetizations] = ...
            SampleGrid(grid,kT,J,numTimePoints,everyT,samplingMethod,timeLag);

%-------------------------------------------------------------------------------
% Plotting:
%-------------------------------------------------------------------------------

f = figure('color','w');
subplot(1,2,1)
imagesc(finalGrid)
axis('square')
colormap([0,0,0;1,1,1])
subplot(1,2,2)
corrMatrix = CorrelationFun(finalGrid,false);
R = RadialAverage(corrMatrix,N);
plot(R,'k')
xlabel('Separation distance')
ylabel('Spatial correlation')



% %------Plot spin-spin correlations------
% f2 = figure(2);
% f2.Color = 'w';

% % Plot spin-spin correlation function
% subplot(121)
% corrMatrix = CorrelationFun(finalGrid,false);
% [X,Y] = meshgrid(-N/2:(N/2-1),-N/2:(N/2-1));
% surf(X,Y,corrMatrix);
% title('Correlation')
% xlabel('x')
% ylabel('y')
% zlabel('Corr')
% colormap('hot')
% caxis([0,min(0.3,max(corrMatrix(:)))])

% % Plot the radial average of spin-spin correlation function
% subplot(122)
% R = RadialAverage(corrMatrix,N);
% plot(R);
% xlabel('Distance')
% ylabel('Correlation')
% title('Radial average')

% %------Plot energy/magnetization variation across Markov chain evolution------
% if exist('energies','var')
%     f3 = figure(3);
%     f3.Color = 'w';
%     subplot(221)
%     plot(energies)
%     xlabel(sprintf('time (iterations/%u)',everyT));
%     ylabel('energy')
%     subplot(222)
%     histogram(energies)
%     xlabel('energy');
%     ylabel('frequency')
%     subplot(223)
%     plot(magnetizations)
%     xlabel(sprintf('time (iterations/%u)',everyT));
%     ylabel('magnetization')
%     subplot(224)
%     histogram(magnetizations)
%     xlabel('magnetization');
%     ylabel('frequency')
% end
