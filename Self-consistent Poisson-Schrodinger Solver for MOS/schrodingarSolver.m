function [psi, en] = schrodingarSolver(m, phi)

%% DEFINE THE QUANTUM REGION
global hcut;
global e;
global dx;
global M;
global N;

X = M - N + 1;        %MAXIMUM DEPLETION REGION LENGTH
V = phi;
%%% OUR SOLUTION WILL BE CONFINED IN THE REGION DEIFNE FROM 
%%% GRID POINT N TO M WHICH DEFINE THE MAXIMUMDEPLETION REGION  
%%% INSIDE SILICON SUBSTRATE

%% FORM THE HAMILTONIAN MATRIX

%%% HERE DIRICHLET BOUNDARY CONDITION IS APPLIED 
%%% THAT IS PSI(0) = 0 AND PSI(N + 1) = 0
H = (diag(ones(1, X) * 2 + V) - diag(ones(1, X - 1), 1) - diag(ones(1, X - 1), -1));
const = ((hcut ^ 2) / (2 * m * dx ^ 2)) / e;
%% SOLVING USING EIGEN SOLVER

[psi, E] = eigs(H, 1, 'SM');
en = E * const;

%%% NORMALIZATION OF WAVE FUNCTION
psi  = psi * sqrt(1 / sum(psi .^ 2));