close all;
clear all;
clc;

%% SEIF-CONSISTANT POISSON-SCHRODINGER SOLVER FOR MOS STRUCTURE
%%% OXIDE - SILICON DI OXIDE; THICKNESS - 50 NM 
%%% SUBSTRATE - P TYPE (BORON DOPED) SILICON; Na - 3 * 10 ^ 16 CM^-3
%%% GATE - ALUMINIUM
%%% DEPLETION TYPE
%%% THRESHOLD VOLTAGE CONSIDERED - 0.65 V
%%% MAXIMUM SPACE CHARGE REGION ~ 200 NM

%%% THE SIMULATION IS DONE OVER THE DOMAIN FROM THE GATE-OXIDE EDGE
%%% TO THE END OF MAXIMUM SPACE CHARGE REGION, THAT IS, A 70 NM RANGE.

%%% AT GATE-OXIDE A DIRICHLET BOUNDARY CONDITION IS USED
%%% AT AT THE END OF SPACE CHARGE REGION A NEUMANN BOUNDARY CONDTION IS USED
global Na;
global tox;
global epSi;
global epOx;
global ep0;

global Nv;
global T;
global e;
global hcut;
global kbT;

global Vg;

Na = (3 * 10 ^ 16) * 10 ^ 6;        % IN M ^ -3
tox = 50 * 10 ^ -9;                 % IN M ^ -3             
epSi = 11.7;
epOx = 3.9;
ep0 = 8.854 * 10 ^ -12;             % IN M ^ -3

Nv = 1.04 * (10 ^ 19) * (10 ^ 6);
T = 300;
e = 1.6 * (10 ^ -19);
hcut = 6.63 * (10 ^ -34) / (2 * pi);
kbT = 0.0259 * (T / 300);

Vg = 0;
%% CALCULATE THE MAXIMUM SPACE CHARGE REGION
global xdT;
global phifp;
ni = (1.5 * 10 ^ 10) * 10 ^ 6;

phifp = 0.0259 * log(Na / ni);
xdT = sqrt(4 * epSi * ep0 * phifp / (e * Na));

%% CALCULATE EFFECTIVE MASS
m0 = 9.1 * (10 ^ -31);

%%% EFFECTIVE MASS FOR CONDUCTIVITY
me = 3 / (1 / 0.89 + 2 / 0.19) * m0;

%%% EFFECTIVE MASS FOR DOS
meDOS = ((6 ^ 2) * 0.89 * 0.19 * 0.19) ^ (1 / 3) * m0;
mhh = 0.49 * m0;
mlh = 0.16 * m0;

%% DEFINING THE MESH
global x;
global M;
global N;
global dx;

M = 101;                           % NUMBER OF GRID POINTS
r = ceil((tox  + xdT) / (10 ^ -9));
x = linspace(0, r, M) * 10 ^ -9;  
dx = x(2) - x(1);                   % MESH SIZE 

% DEPLETION REGION START FROM N
for j = 1:M
    disp(x(j) - tox)
    if x(j) >= tox
        N = j;
        break
    end
end

%% INITIALIZER
%%% CREATE A INITIAL POTENTIAL PROFILE TO WORK ON
v = ones(1, M);
v(1:N) = 3.25;
v(N:M) = 0;
phi = v(N:M);

for i=1:300
%% SCHRODINGAR SOLVER
[psi, en] = schrodingarSolver(me, phi); 

%% POISSON SOLVER
[phi] = poissonSolver();
phi = phi + 0.6 * (phinew - phi); 


end

figure(1)
plot(x(N:M), phinew);

figure(2)
plot(x(N:M), psi)