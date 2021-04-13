%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% SCRIPT TO COMPUTE HOW THE SPATIAL CORRELATION LENGTH CHANGES WITH kT
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% SET PARAMETERS
%-------------------------------------------------------------------------------
% N, linear lattice size
N = 100;
% J, coupling strength (change sign for antiferromagnetic coupling!)
J = 1;
% numTimePoints, number of update steps (use large multiple of N^2 for Metropolis)
numTimePoints = 1000*N^2;
% everyT, plot and store the energy/magnetization of the grid everyT iterations
everyT = 100*N^2;
% reInitialize, whether to generate a new initial condition (or continue from previous)
reInitialize = true;
% p, average proportion of initial +1 spins
p = 0.8; % (0.5 for random initial condition)
% samplingMethod, 'HeatBath', 'Metropolis' or 'Wolff'
samplingMethod = 'Metropolis';
% timeLag
timeLag = 0; % option to slow down plotting

%-------------------------------------------------------------------------------
% Run the sampling algorithm
%-------------------------------------------------------------------------------
kT = 2:0.02:2.5;
numkT = length(kT);
numRepeats = 5;
lambda = zeros(numkT,numRepeats);
for i = 1:numkT
    for r = 1:numRepeats
        grid = sign(p-rand(N)); % random initial configuration
        [finalGrid,energies,magnetizations] = ...
                SampleGrid(grid,kT(i),J,numTimePoints,everyT,samplingMethod,timeLag);
        corrMatrix = CorrelationFun(finalGrid,false);
        R = RadialAverage(corrMatrix,N);
        firstDrop = find(R < 1/exp(1),1,'first');
        if isempty(firstDrop)
            lambda(i,r) = length(R);
        else
            lambda(i,r) = firstDrop;
        end
        fprintf(1,'kT = %.2f, lambda = %.0f\n',kT(i),lambda(i));
    end
end

%-------------------------------------------------------------------------------
lambdaMean = mean(lambda,2);
lambdaStd = std(lambda,0,2);
%-------------------------------------------------------------------------------
f = figure('color','w');
hold('on')
% plot(kT,lambdaMean,'o-k')
plot([min(kT),max(kT)],ones(2,1)*N/2,'--','color',ones(3,1)*0.5,'LineWidth',1)
plot(2/log(1+sqrt(2))*ones(2,1),[0,N/2],'LineWidth',1)
errorbar(kT,lambdaMean,lambdaStd,'o-k','LineWidth',1.5)
xlabel('kT')
ylabel('Correlation length')
title(sprintf('%u x %u lattice',N,N))
ylim([0,N/2])
