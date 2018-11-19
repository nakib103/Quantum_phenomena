close all;
clear all;
clc;

%% SEIF-CONSISTANT POISSON-SCHRODINGER SOLVER FOR MOS STRUCTURE
%%% OXIDE: SILICON DI OXIDE; THICKNESS - 1 NM 
%%% SUBSTRATE: P TYPE SILICON; Na - 10 ^ 18 CM^-3
%%% GATE: ALUMINIUM
%%% DEPLETION TYPE

Na = (10 ^ 18) * (10 ^ 6);          % doping level
Nv = 1.04 * (10 ^ 19) * (10 ^ 6);
ep0 = 8.854 * (10 ^ -12);            
epSi = 11.7 * ep0;                  % permittivity of silicon
epOx = 3.9 * ep0;                   % permittivity of oxide

T = 300;                            
e = 1.6 * (10 ^ -19);               
hcut = 6.63 * (10 ^ -34) / (2 * pi);
kbT = 0.0259 * (T / 300);

m0 = 9.1 * (10 ^ -31);
mC = 1.08 * m0;                     % conduction bend effective mass
mDOS = 0.26 * m0;                   % density of states effective mass

%% DEFINING THE MESH
nm = 10 ^ -9;
dx = 0.1 * nm;                      % mesh size

tox = (1) * nm;                     % oxide thickness
tsub = (50) * nm;                   % substrate thickness
x = -20*nm:dx:tsub;
len = length(x);

vg = 0;
phim = 3.2;                         % modified metal work function
phiX = 3.25;                        % modified oxide electron affinity

%% CALCULATE THE MAXIMUM SPACE CHARGE REGION
ni = (1.5 * 10 ^ 10) * 10 ^ 6;
phifp = 0.0259 * log(Na / ni);

Vox0 = 0.05;                        
surPot= abs(3.2 - (3.25 + 0.56 + phifp) - Vox0);    % surface potantial
xd = sqrt((2 * epSi * surPot) / (e * Na));          % deplation region

vf = (0.56 + phifp);

%% POISSON EQUATION ON OXIDE
for i = 1:len;
    if x(i) >= -tox;
        oxi = i;
        break
    end
end
for i = oxi:len;
    if (x(i) >= 0)
        oxj = i;
        break;
    end
end
LAPox = (-diag(ones(1,oxj-oxi+1) * 2)) + (diag(ones(1,oxj-oxi),1)) + (diag(ones(1,oxj-oxi),-1));

B = zeros(oxj-oxi+1, 1);
B(1) = -phim;
B(end) = - (abs(vf - surPot) + phiX);

vOx = inv(LAPox) * B;

%% POISSON EQUATION ON INVERSION REGION
for i = 1:len;
    if x(i) >= xd;
        depj = i;
        break;
    end
end
depi = oxj + 1;

LAPdep = (-diag(ones(1,depj-depi+1) * 2)) + (diag(ones(1,depj-depi),1)) + (diag(ones(1,depj-depi),-1));

B = -ones(depj-depi+1, 1)*(e * Na * dx^2) / epSi;
B(1) = B(1) - (vOx(end) - phiX);
B(end) = B(end) - vf;

vdep = inv(LAPdep) * B;

%% INITIAL TRIAL VOLTAGE
vmetal = zeros(oxi-1,1);
vsub = ones(len-depj,1) * vf;
Vtot = [vmetal' vOx' vdep' vsub'];

figure(1)
plot(x/nm, Vtot)
xlabel('device length (nm)');
ylabel('trial potential (eV)');

%% ITERATION
const = ((hcut ^ 2) / (2 * mC * dx ^ 2)) / e;
lendep = depj - depi + 2;
V = [vOx(end) vdep'];
enBound = zeros(1, lendep);
psiBound = zeros(lendep, lendep);
Ec_Ef = 1.12 - kbT * log(Nv / Na);
n = zeros(lendep, 1);

error = zeros(1, 100);
for iter = 1:100

%% SCHRODINGAR SOLVER
%%% FORM THE HAMILTONIAN MATRIX
DEL = (-diag(ones(1,lendep) * 2)) + (diag(ones(1,lendep-1),1)) + (diag(ones(1,lendep-1),-1));
H = -(DEL * const) + (diag(V'));
%%% SOLVING USING EIGEN SOLVER
[psi, en] = eig(H);

en = diag(en);
bound = 1;

for j = 1:lendep
    if (en(j) > vdep(1) && en(j) < vf)
        enBound(bound) = en(j);
        psiBound(:, bound) = psi(:, j);
        bound = bound + 1;
    end
end

%%% CALCULATING CARRIER DISTRIBUTION
Ef_Ei = -(Ec_Ef - surPot + enBound);     
Ni = (e * kbT * mDOS) / (pi * hcut^2) * log(1 + exp(Ef_Ei/kbT));
psi2 = (psiBound.*psiBound);    
k = Ni ./ sum(psi2);

for j = 1:bound-1
    n = n + k(j) * psi2(:,j);
end

%% POISSON SOLVER
LAPdep = (-diag(ones(1,lendep) * 2)) + (diag(ones(1,lendep-1),1)) + (diag(ones(1,lendep-1),-1));
B = - (Na - n) * (e * dx^2 / epSi);
B(1) = B(1) - (vOx(end) - phiX);
B(end) = B(end) - surPot;

Vnew = inv(LAPdep) * B;

err = sqrt(sum((Vnew' - V).^ 2) / lendep);
disp(iter)
disp(err)
error(iter) = err;

V = V + 0.3 * (Vnew' - V);

if (err < 10^-9)
    break
end

if iter == 1
    figure(2)
    plot((1:lendep)*dx/nm, psi(:,1).*psi(:,1), 'r')
    hold on
    plot((1:lendep)*dx/nm, psi(:,2).*psi(:,2), 'b')
    hold on
    plot((1:lendep)*dx/nm, psi(:,3).*psi(:,3), 'g')
    
    xlabel('device length (nm)');
    ylabel('wave function');
end

if iter == 50
    figure(4)
    plot((1:lendep)*dx/nm, psi(:,1).*psi(:,1), 'r')
    hold on
    plot((1:lendep)*dx/nm, psi(:,2).*psi(:,2), 'b')
    hold on
    plot((1:lendep)*dx/nm, psi(:,3).*psi(:,3), 'g')
    
    xlabel('device length (nm)');
    ylabel('wave function');
end

figure(5)
plot((1:lendep)*dx/nm, V)
hold on
xlabel('device length (nm)');
ylabel('potential (eV)');

figure(6)
plot((1:lendep)*dx/nm, n)
hold on

xlabel('device length (nm)');
ylabel('bound carrier distribution (m^-2)');
end

vsub = ones(1, length(vsub))*V(end);
Vtot = [vmetal' vOx' V(2:end) vsub];
figure(7)
plot(x/nm, Vtot)
xlabel('device length (nm)');
ylabel('potential (eV)');