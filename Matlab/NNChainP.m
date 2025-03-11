function h = NNChainP(a,b,varargin)
% Defines a Hamiltonian with a(k)b(k+1) nearest-neighbor
% and periodic boundary condition.
% NNChainP(a,b,n) defines a a(k)b(k+1) type Hamiltonian
% with periodic boundary conditions.
% If argument n is omitted than the default is taking to be
% the value of global variable N.
%-------------------------------------------------------------------------------

if isempty(varargin),
    global N;
else
    if length(varargin)~=1,
        error('Wrong number of input arguments.')
    end %if
    N = varargin{1};
end

h = sparse(2^N,2^N);
op = kron(a,b);
for n = 1:N-1
    h = h+kron(kron(speye(2^(n-1)),op),speye(2^(N-n-1)));
end

% Periodic boundary conditions
if N > 1
   h = h+kron(kron(b,speye(2^(n-1))),kron(speye(2^(N-n-1)),a));
end
