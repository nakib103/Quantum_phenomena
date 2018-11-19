close all; 
clear all;
clc;

% constants
hcut =  6.63 * (10 ^ -34) / (2 * pi);
m = 9.1 * (10 ^ -31);
e = 1.6 * (10^ -19);
e0 = 8.854 * (10 ^ -12);


% defining the mesh
nm = 10^-9;
dx = 0.00625 * nm;                  % mesh size
x1 = -0.125*nm:dx:0.125*nm-dx;         % defining the x1 dimension
N = length(x1);
    

k = - hcut^2 / (2 * m * e * dx^2);
C = - e / (4 * pi * e0);
%% defining the kinetic energy matrix
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
%% defining the potential well
A = ones(1,N);
for i = 2:N
    A(i) = 1 / ((i-1)*dx);
end
A(1) = 1 / 10^-20;

V = zeros(N^2, N^2);
for i = 1:N
    K = circshift(A, [0, (i - 1)]);
    for j = i-1:-1:1
        K(j) = 1 / ((i-j)*dx);
    end
    K = diag(K);
    V((i-1)*N + 1:(i)*N, (i-1)*N + 1:(i)*N) = K;
end

V = C*V;
%% eigensolver 
H = H + V;
[WF, En] = eig(H);         

%% GRAPHS
figure(1)
psi = reshape(WF(:, N^2), N, N);

surf(psi.^2)
%view([0 90])
xlabel('x1')
ylabel('x2')
