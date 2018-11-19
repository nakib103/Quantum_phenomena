function [phi] = poissonSolver(me, mhh, mlh, en, psi)

%% FORMING THE RHS MATRIX
global M;
global N;
global e;
global epSi;
global epOx;
global ep0;
global dx;
global hcut;
global Vg;

%%% CALCULATE THE IONIZED CHARGE CARRIERS
global Na;          % DONOR ATOMS
global kbT;           
global Nv;
G = 4;              % FOR HEAVY AND LIGHT HOLE DEGENERACY IS G = 4
ioen = 0.045;       % IONIZATION ENERGY OF BORON IS 45 MEV                
ionizedNa = Na / (1 + G * exp((0.506 - ioen) / kbT));

%%% CALCULATE THE MOBILE CHARGE CARRIERS
%%% CALCULATE FERMI LEVEL POSITION
Ef = kbT * log(Nv / Na);

%%% CALCULATE EFFECTIVE DOS IN 2D
dose = (me * kbT * e) / (pi * hcut ^ 2);
%doshh = (mhh * kbT) / (pi * hcut ^ 2);
%doslh = (mlh * kbT) / (pi * hcut ^ 2);

%%% CALCULATE FERMI INTEGRAL
etaF = (Ef - en) / kbT;
fermiInge = log(1 + etaF);

%etaF = (0 - Ef) / kbT;
%fermiIngh = log(1 + etaF);

%%% CALCULATE DENSITY
n = dose * fermiInge * (psi.*psi);
%p = (doshh + doslh) * fermiIngh;

%%% TOTAL DENSITY
rho = ionizedNa - n';

%%% RHS MATRIX
%%% NO CHARGE CONSIDERED IN OXIDE (i.e. SURFACE CHARGE NOT CONSIDERED)
B = ones(1, M-N+1).* ((rho * e) / (epSi * ep0)) ;
size(B)
global phifp
vf = (phifp + 3.25) - 3.2;
B(end) = B(end) - vf;
%% FORMING THE LHS MATRICES

%ep = ones(1, M + 1);
%ep(1:N-1) = ep0 * epOx;
%ep(N:M+1) = ep0 * epSi;

%A = ( -(diag(ones(1,M).* (ep(1:M) + ep(2:M+1)) )) + (diag(ones(1,M-1) .* ep(2:M),1))...
%    + (diag(ones(1,M-1) .* ep(1:M-1),-1))) / dx^2;
A = ( -(diag(ones(1,M-N+1) * 2)) + (diag(ones(1,M-N),1)) + (diag(ones(1,M-N),-1)) ) / dx^2;

%%% ADDING ANOTHER MATRIX FOR BOUNDARY CONDITION
%%% WE USE DIRICHLET BOUNDARY CONDITION

%%% AT GATE-SIO2 INTERFACE PHI = GATE VOLTAGE - WORK
%%% FUNCTION DIFFERENCE OF METAL AND SIO2

%%% AT THE EDGE OF THE SPSCE CHARGE REGION WE SET PHI = 0
%D = zeros(1, M);
%D(1) = -(Ef - Vg - 3.2) * ep(1) / dx^2;
%D(M) = -1.12 * ep(M+1) / dx^2;
size(B)
%% SOLUTION
C = inv(A)*(B)';
phi = C';
figure(3)
plot(phi)