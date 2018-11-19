function [phi] = initialize()

global N;
global M;

phi = ones(1, M);

x = N:M;
y = 1 ./ (1 + exp(-x));

phi(N:M) = y;