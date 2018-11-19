close all; 
clear all;
clc;

% constants
hcut =  6.63 * (10 ^ -34) / (2 * pi);
m = 9.1 * (10 ^ -31);
e = 1.6 * (10^-19);

% defining the mesh
nm = 10^-9;                       
dx = 0.1 * nm;                  % mesh size
x1 = -2*nm:dx:2*nm-dx;             % defining the x1 dimension (-2nm - 2nm)
N = length(x1);
    

k = - hcut^2 / (2 * m * e * dx^2);
%% defining the kinetic energy matrix
%%% used dirichlet boundary condition
B = (diag(ones(1,N)*4)-diag(ones(1,N-1),1)-diag(ones(1,N-1),-1));
I = -eye(N, N);
H = zeros(N^2, N^2);

for i = 2:N-1
    H((i-1)*N + 1:(i)*N, (i-2)*N + 1:(i-1)*N) = I;
    H((i-1)*N + 1:(i)*N, (i-1)*N + 1:(i)*N) = B;
    H((i-1)*N + 1:(i)*N, (i)*N + 1:(i+1)*N) = I;
end
H(1:N,1:N) = B;
H(1:N, N+1:2*N) = I;
H((N-1)*N + 1:N^2,(N-2)*N + 1:(N-1)*N) = I;
H((N-1)*N + 1:N^2,(N-1)*N + 1:N^2) = B;

H = k*H;

%% eigensolver 
[WF, En] = eig(H);         

%% GRAPHS
psi = reshape(WF(:, N^2), N, N);
psi1 = reshape(WF(:, N^2-1), N, N);
psi2 = reshape(WF(:, N^2-2), N, N);
psi3 = reshape(WF(:, N^2-3), N, N);
psi4 = reshape(WF(:, N^2-4), N, N);
psi5 = reshape(WF(:, N^2-5), N, N);
psi6 = reshape(WF(:, N^2-6), N, N);
psi7 = reshape(WF(:, N^2-7), N, N);

figure(1)
surf(psi.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')

figure(2)
subplot(121)
surf(psi1.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')

subplot(122)
surf(psi2.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')


figure(4)
surf(psi3.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')

figure(5)
subplot(121)
surf(psi4.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')

subplot(122)
surf(psi5.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')

figure(6)
subplot(121)
surf(psi6.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')

subplot(122)
surf(psi7.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')
%% energy
E = diag(En);