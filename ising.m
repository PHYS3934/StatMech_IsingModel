% Choose kT, N (linear lattice size), and J (strength of the coupling)
kT = 2/log(1+sqrt(2));
N = 50;
J = 1; % change sign for antiferromagnetic coupling

% Generate a random initial configuration. 
% Comment out to keep sampling with the previous configuration.
p=.5; % average proportion of initial +1 spins
grid = sign(p-rand(N)); % random initial configuration

% Run the Wolff or Metropolis algorithm and return a matrix of spin values.
t =10*N^2;% choose t update steps (use big multiple of N^2 for Metropolis)
grid = metropolis(N,kT,J,t,grid);
%grid = wolff(N,kT,J,t,grid);

% Compute final magnetization density and energy density
M = sum(sum(grid))/numel(grid);
E = isingenergy(grid,J);

% Plot correlation function
cor = correlation(grid); figure(2); surf(cor); 
R = radialavg(cor,N); figure(3); plot(R);